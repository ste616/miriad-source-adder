#!/usr/bin/python
"""Miriad Source Adding Helper

Usage:
  miriad-source-adder.py [--source-file=<str>] [--ra=<str> ...] [--dec=<str> ...] [--flux=<flux> ...] [--size=<bmaj,bmin,bpa> ...] [--alpha=<alpha> ...] [--out=<str>] [--test] [--temp-dir=<str>] <dataset>

-h --help                show this
-s --source-file FILE    file to give to uvgen as "source" parameter
-r --ra RA               the RA of a source to add HH:MM:SS.S
-d --dec DEC             the declination of a source to add DDD:MM:SS.S
-f --flux FLUX           the flux density of a source to add in Jy
-S --size BMAJ,BMIN,BPA  the size of a source to add, in arcsec for each length, and degrees for the angle
-a --alpha ALPHA         the spectral index of a source to add
-o --out OUT             output the mixed dataset with this name
-t --test                only do a single source uvgen for testing purposes
-T --temp-dir DIR        put all the intermediate files in this directory [default: .]
"""

from docopt import docopt
import math
from mirpy import miriad
import ephem
from datetime import date, datetime
import os
import shutil
import numpy as np
import sys

# A routine to turn a Miriad type time string into a ephem Date.
def mirtime_to_date(mt):
    year = 2000 + int(mt[0:2])
    monthShort = mt[2:5]
    date = int(mt[5:7])
    hour = int(mt[8:10])
    minute = int(mt[11:13])
    second = int(round(float(mt[14:18])))
    monthDict = { 'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
                  'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12 }
    dateString = "%4d/%02d/%02d %02d:%02d:%02d" % (year, monthDict[monthShort], date,
                                                   hour, minute, second)
    return ephem.Date(dateString)

def round_date_5second(dt):
    # This routine goes to the nearest time that ends in a 5 second.
    d = dt.datetime()
    dsecd = d.second % 10
    dd = 5 - dsecd
    dt = ephem.date(dt + dd * ephem.second)
    return dt

def round_date_nsecond(dt, n):
    # This routine takes some date, and returns the date nearest to
    # that which falls an integer number of n-second intervals away
    # from midnight on that day.
    d = dt.datetime()
    # How many seconds after midnight do we have?
    ts = d.second + d.minute * 60 + d.hour * 3600
    # Work out the nearest cycle centre.
    ccs = round(((ts - n / 2) / n)) * n + (n / 2)
    dd = ccs - ts
    rt = ephem.date(dt + dd * ephem.second)
    return rt

def date_to_mirtime(dt):
    # Output a Miriad formatted date.
    d = dt.datetime()
    rs = d.strftime("%y%b%d:%H:%M:%S").lower()
    return rs


# Use uvindex to work out the necessary parameters of this dataset.
def filter_uvindex(output):
    # We send back a dictionary.
    rd = { 'index': { 'time': [], 'source': [], 'calcode': [], 'antennas': [],
                      'spectral_channels': [], 'wideband_channels': [], 'freq_config': [],
                      'record_number': [] },
           'total_time': 0,
           'freq_configs': [], 'polarisations': [], 'sources': [] }
    outlines = output.split('\n')
    section = 0
    freqconfig_n = 0
    freqconfig_found = 0
    fc = None
    sourcearea = 0
    for i in xrange(0, len(outlines)):
        index_elements = outlines[i].split()
        if ((section == 0) and (len(outlines[i]) >= 74)):
            indexTime = mirtime_to_date(outlines[i][0:18])
            if ((index_elements[1] != "Total") and (index_elements[2] != "number")):
                # This is a regular line.
                offset = 0
                rd['index']['time'].append(indexTime)
                rd['index']['source'].append(index_elements[1])
                # Check if we have a calibrator code.
                calcode = outlines[i][36:37]
                if (calcode == " "):
                    # No calibrator code.
                    offset = 1
                rd['index']['calcode'].append(calcode)
                rd['index']['antennas'].append(int(index_elements[3 - offset]))
                rd['index']['spectral_channels'].append(int(index_elements[4 - offset]))
                rd['index']['wideband_channels'].append(int(index_elements[5 - offset]))
                rd['index']['freq_config'].append(int(index_elements[6 - offset]))
                rd['index']['record_number'].append(int(index_elements[7 - offset]))
            else:
                # We've moved to the next section
                section = 1
        elif ((section == 1) and (len(index_elements) > 0) and (index_elements[0] == "Total") and
              (index_elements[1] == "observing")):
            # We've found the total amount of observing time.
            rd['total_time'] = float(index_elements[4])
            section = 2
        elif (section == 2):
            if ((len(index_elements) > 0) and (index_elements[0] == "Frequency")
                and (index_elements[1] == "Configuration")):
                freqconfig_n = int(index_elements[2])
                freqconfig_found = 1
                if (fc is not None):
                    rd['freq_configs'].append(fc)
                fc = { 'number': freqconfig_n, 'nchannels': [],
                       'frequency1': [], 'frequency_increment': [],
                       'rest_frequency': [], 'ifchain': [] }
            elif (freqconfig_found == 1):
                freqconfig_found = 2
            elif (freqconfig_found == 2):
                if (outlines[i] == ""):
                    freqconfig_found = 0
                else:
                    # This is the actual line.
                    fc['nchannels'].append(int(index_elements[0]))
                    fc['frequency1'].append(float(index_elements[1]))
                    fc['frequency_increment'].append(float(index_elements[2]))
                    fc['rest_frequency'].append(float(index_elements[3]))
                    fc['ifchain'].append(int(index_elements[5]))
            elif (outlines[i] == "------------------------------------------------"):
                if (fc is not None):
                    rd['freq_configs'].append(fc)
                section = 3
        elif (section == 3):
            if ((len(index_elements) > 0) and (index_elements[0] == "There") and
                (index_elements[3] == "records") and (index_elements[5] == "polarization")):
                rd['polarisations'].append(index_elements[6])
            elif (outlines[i] == "------------------------------------------------"):
                section = 4
        elif (section == 4):
            if ((len(index_elements) > 0) and (index_elements[0] == "Source")):
                sourcearea = 1
            elif ((len(index_elements) > 2) and (sourcearea == 1)):
                src = { 'name': index_elements[0], 'calcode': index_elements[1],
                        'right_ascension': index_elements[2], 'declination': index_elements[3],
                        'dra': index_elements[4], 'ddec': index_elements[5] }
                rd['sources'].append(src)

    # Convert things into numpy arrays for easy where-ing later.
    rd['index']['time'] = np.array(rd['index']['time'])
    rd['index']['source'] = np.array(rd['index']['source'])
    rd['index']['calcode'] = np.array(rd['index']['calcode'])
    rd['index']['antennas'] = np.array(rd['index']['antennas'])
    rd['index']['spectral_channels'] = np.array(rd['index']['spectral_channels'])
    rd['index']['wideband_channels'] = np.array(rd['index']['wideband_channels'])
    rd['index']['freq_config'] = np.array(rd['index']['freq_config'])
    rd['index']['record_number'] = np.array(rd['index']['record_number'])
    
    return rd

# Use uvlist to get the actual positions of the antennas.
def filter_uvlist_antennas(output):
    # We send back a dictionary.
    rd = { 'telescope': "", 'latitude': "", 'longitude': "", 'antennas': [] }
    outlines = output.split('\n')
    coords = 0
    for i in xrange(0, len(outlines)):
        index_elements = outlines[i].split()
        if (len(index_elements) > 0):
            if (index_elements[0] == "Telescope:"):
                rd['telescope'] = index_elements[1]
            elif (index_elements[0] == "Latitude:"):
                rd['latitude'] = index_elements[1]
            elif (index_elements[0] == "Longitude:"):
                rd['longitude'] = index_elements[1]
            elif ((len(index_elements) == 3) and (index_elements[1] == "----------")):
                coords = 1
            elif (coords == 1):
                ant = { 'number': int(index_elements[0]), 'coord_x': float(index_elements[1]),
                        'coord_y': float(index_elements[2]), 'coord_z': float(index_elements[3]) }
                rd['antennas'].append(ant)
    return rd

# Use a uvlist log to get the cycle time.
def filter_uvlist_variables(logfile_name):
    # We send back a dictionary.
    rd = { 'cycle_time': -1. }
    with open(logfile_name, "r") as fp:
        loglines = fp.readlines()
    for i in xrange(0, len(loglines)):
        index_elements = loglines[i].split()
        if ((len(index_elements) > 2) and
            (index_elements[0] == "inttime") and (index_elements[1] == ":")):
            rd['cycle_time'] = float(index_elements[2])
    return rd

def get_hour_angle(source, observer):
    source.compute(observer)
    hour_angle = (observer.sidereal_time() - source._ra) * 180. / (math.pi * 15.)
    if (hour_angle > 12):
        hour_angle -= 24.
    if (hour_angle < -12):
        hour_angle += 24.
    return hour_angle

def output_antenna_file(location, filename):
    lines = []
    for i in xrange(0, len(location['antennas'])):
        lines.append("%12.4f %12.4f %12.4f" % (location['antennas'][i]['coord_x'],
                                               location['antennas'][i]['coord_y'],
                                               location['antennas'][i]['coord_z']))
    with open(filename, "w") as fp:
        for i in xrange(0, len(lines)):
            fp.write("%s\n" % lines[i])

def split_into_segments(idx):
    # We go through a uvindex dictionary and return segments.
    # Each segment is a single source, with a start and end time.
    segs = []
    oldsrc = ""
    sseg = None
    for i in xrange(0, len(idx['index']['source'])):
        if (idx['index']['source'][i] != oldsrc):
            if (oldsrc != ""):
                # Put the segment on the list.
                segs.append(sseg)
            oldsrc = idx['index']['source'][i]
            sseg = { 'source': idx['index']['source'][i],
                     'start_time': ephem.Date(idx['index']['time'][i]),
                     'end_time': ephem.Date(idx['index']['time'][i]) }
        else:
            sseg['end_time'] = idx['index']['time'][i]
    # Have to push the last segment on.
    segs.append(sseg)
    return segs

def add_uvcat(fileString, newFile, outdir):
    # Keep adding files to the uvcat input string until it gets too long
    # and then do the uvcat.
    lfs = len(fileString)
    lnf = len(newFile)
    maxlen = 950
    if ((lfs + lnf + 1) < maxlen):
        # Easy, just append the string.
        return "%s,%s" % (fileString, newFile)
    # Otherwise we do the uvcat.
    files = fileString.split(",")
    ohead, otail = os.path.split(files[0])
    outFile = "%s/uvc_%s" % (outdir, otail)
    if (os.path.isdir(outFile)):
        shutil.rmtree(outFile)
    print "  Concatenating intermediate product"
    #print "[%s]" % fileString
    miriad.uvcat(vis=fileString, out=outFile)
    # Now add the new file to this.
    rv = add_uvcat(outFile, newFile, outdir)
    return rv

def stringToFloat(s):
    # Take a sexagesimal string and return the float value.
    cmps = s.split(":")
    degs = float(cmps[0])
    neg = False
    if (cmps[0][0] == "-"):
        neg = True
    if neg:
        degs = degs * -1.
    mins = float(cmps[1])
    secs = float(cmps[2])
    val = degs + mins / 60. + secs / 3600.
    if neg:
        val = val * -1.
    return val

def add_source(args):
    # Let's create a uvindex of the dataset first.
    dataset_name = args['<dataset>']
    miriad.set_filter('uvindex', filter_uvindex)
    index_data = miriad.uvindex(vis=dataset_name, interval="0.1")

    # Get the cycle time.
    uvlist_log_name = "uvlist.log"
    if (os.path.isfile(uvlist_log_name)):
        os.remove(uvlist_log_name)
    # We write to a log because otherwise this command waits for the user to press
    # the enter key halfway through operation.
    miriad.uvlist(vis=dataset_name, options="variables,full",
                  log=uvlist_log_name)
    telescope_variables = filter_uvlist_variables(uvlist_log_name)
    cycle_time = round(telescope_variables['cycle_time'])
    print "     Found cycle time of %d seconds" % cycle_time
    
    # Get the position of the telescope.
    miriad.set_filter('uvlist', filter_uvlist_antennas)
    telescope_coordinates = miriad.uvlist(vis=dataset_name, options="array,full")
    
    # Set up a pyephem observatory for this telescope.
    telescope = ephem.Observer()
    telescope.pressure = 0
    telescope.horizon = '12'
    telescope.lat = telescope_coordinates['latitude']
    telescope.lon = telescope_coordinates['longitude']

    # We generate a dataset for each source that we found.
    # In general, if there is more than one source in the dataset, it's probably
    # a mosaic.

    # First, print a warning if the user has specified a source file (which only
    # supports offsets) if there is more than one source.
    if ((len(index_data['sources']) > 1) and (arguments['--source-file'] is not None)):
        print "WARNING: multiple sources present, but a source file has been given."
        print "         This is probably not what you wanted, but we continue in case it is."
    

    # We'll go through in segments, as this will handle mosaic observations
    # more efficiently.
    segments = split_into_segments(index_data)

    # Work out how many of these things to do.
    num_gens = len(segments)
    if (args['--test']):
        num_gens = 1

    # The string to keep as input for the uvcat-ing.
    uvcatString = ""
    for i in xrange(0, num_gens):
        # Make a source object for this position.
        source = ephem.FixedBody()
        source.name = segments[i]['source']
        source._epoch = "2000"
        source_index = -1
        for j in xrange(0, len(index_data['sources'])):
            if (index_data['sources'][j]['name'] == segments[i]['source']):
                source._ra = index_data['sources'][j]['right_ascension']
                source._dec = index_data['sources'][j]['declination']
                source_index = j
                break
        source.compute()
        sourceRa = 15. * stringToFloat(index_data['sources'][j]['right_ascension'])
        sourceDec = stringToFloat(index_data['sources'][j]['declination'])
        
        
        ##### Find the times that this source was observed.
        # Where was the source at the start of the observations?
        # Go to the starting time of the observation.
        telescope.date = ephem.Date(segments[i]['start_time'])
        # We really just need the hour angle.
        start_hour_angle = get_hour_angle(source, telescope)

        # And the same for the end of the observation.
        telescope.date = ephem.Date(segments[i]['end_time'])
        finish_hour_angle = get_hour_angle(source, telescope)
        #print "%s %f %f" % (source.name, start_hour_angle, finish_hour_angle)

        # Work out the 0 hour angle. Get the nearest time to 0.
        transit_time = None
        if (abs(start_hour_angle) < abs(finish_hour_angle)):
            transit_time = segments[i]['start_time'] - (start_hour_angle * ephem.hour)
        else:
            transit_time = segments[i]['end_time'] - (finish_hour_angle * ephem.hour)
        transit_time = ephem.Date(transit_time)

        #print "transit time is %s" % transit_time

        # Work out the inputs to uvgen.
        if (arguments['--source-file'] is not None):
            uvgen_source = arguments['--source-file']
        else:
            # We create the source file per pointing.
            # Determine the offset between the sources we want to add and the pointing centre.
            src_lines = []
            for j in xrange(0, len(args['--ra'])):
                nsource = ephem.FixedBody()
                nsource._epoch = "2000"
                nsource._ra = args['--ra'][j]
                nsource._dec = args['--dec'][j]
                nsource.compute()
                nsourceRa = 15. * stringToFloat(args['--ra'][j])
                nsourceDec = stringToFloat(args['--dec'][j])
                #dra = (nsource.ra - source.ra) * (180. * 3600.) / math.pi # in arcseconds
                #ddec = (nsource.dec - source.dec) * (180. * 3600.) / math.pi # also in arcseconds
                avdec = ((nsourceDec + sourceDec) / 2.) * (math.pi / 180.)
                dra = (nsourceRa - sourceRa) * 3600. * math.cos(avdec)
                ddec = (nsourceDec - sourceDec) * 3600.
                # The source flux density.
                fluxdens = float(args['--flux'][j]) # in Jy
                # The source size.
                szel = args['--size'][j].split(",")
                if (len(szel) != 3):
                    print "ERROR: the size of source %d was not correctly specified" % (j + 1)
                    sys.exit()
                bmaj = float(szel[0])
                bmin = float(szel[1])
                bpa = float(szel[2])
                # The spectral index.
                specidx = float(args['--alpha'][j])
                src_lines.append("%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,0,0,0,%.3f" % (fluxdens, dra, ddec,
                                                                               bmaj, bmin, bpa, specidx))
            # Make the file.
            uvgen_source = "%s/source_created_%s" % (args['--temp-dir'],
                                                     source.name)
            #print "  Creating source generation file %s" % uvgen_source
            with open(uvgen_source, "w") as fp:
                for j in xrange(0, len(src_lines)):
                    fp.write("%s\n" % src_lines[j])
            
        uvgen_telescop = telescope_coordinates['telescope'].lower()
        uvgen_stokes = ",".join(index_data['polarisations']).lower()
        #print uvgen_source
        #print uvgen_telescop
        #print uvgen_stokes
        uvgen_lat = "%.6f" % math.degrees(telescope.lat)
        #print uvgen_lat
        uvgen_radec = "%s,%s" % (index_data['sources'][source_index]['right_ascension'],
                                 index_data['sources'][source_index]['declination'])
        #print uvgen_radec
        # Extend the HA range a little on each side.
        sidereal_modifier = 1.00273790935
        start_ha = math.floor(start_hour_angle * 10. * sidereal_modifier) / 10. - 0.05
        finish_ha = math.ceil(finish_hour_angle * 10. * sidereal_modifier) / 10. + 0.05
        uvgen_harange = "%.2f,%.2f" % (start_ha, finish_ha)
        #print uvgen_harange
        # The time at 0 HA.
        uvgen_time = date_to_mirtime(transit_time)
        #print uvgen_time
        # Use only frequency config 0, but we could probably change that if we need to.
        #print index_data['freq_configs']
        chan_offset = math.floor(index_data['freq_configs'][0]['nchannels'][0] / 2)
        chan_spacing_mhz = index_data['freq_configs'][0]['frequency_increment'][0] * 1000.
        freq_offset_mhz = chan_spacing_mhz * chan_offset
        width_mhz = chan_spacing_mhz * index_data['freq_configs'][0]['nchannels'][0]
        uvgen_corr = "%d,1,%d,%.3f" % (index_data['freq_configs'][0]['nchannels'][0],
                                       freq_offset_mhz, width_mhz)
        #print uvgen_corr
        uvgen_freq = "%.3f,0.0" % index_data['freq_configs'][0]['frequency1'][0]
        #print uvgen_freq
        # Now we have to write out the antenna locations file.
        uvgen_ant = "antenna_configuration.file"
        output_antenna_file(telescope_coordinates, uvgen_ant)
        uvgen_baseunit = 3.33564
        
        # The name of the simulated dataset.
        simulated_name = "%s/%s_%s.uvgen" % (args['--temp-dir'], args['--out'],
                                             source.name)
        # Delete this set if it exists.
        if os.path.isdir(simulated_name):
            shutil.rmtree(simulated_name)

        # Run uvgen now.
        print "  Running uvgen for source %d / %d" % ((i + 1), num_gens)
        with open("last_uvgen.dbg", "w") as fp:
            fp.write("uvgen source=%s ant=%s baseunit=%s telescop=%s corr=%s time=%s freq=%s radec=%s harange=%s stokes=%s lat=%s out=%s inttime=%s\n" % (uvgen_source, uvgen_ant, uvgen_baseunit, uvgen_telescop, uvgen_corr, uvgen_time, uvgen_freq, uvgen_radec, uvgen_harange, uvgen_stokes, uvgen_lat, simulated_name, cycle_time))
        miriad.uvgen(source=uvgen_source, ant=uvgen_ant, baseunit=uvgen_baseunit,
                     telescop=uvgen_telescop, corr=uvgen_corr, time=uvgen_time,
                     freq=uvgen_freq, radec=uvgen_radec, harange=uvgen_harange,
                     stokes=uvgen_stokes, lat=uvgen_lat, out=simulated_name,
                     inttime=cycle_time)
        # Chop out just the time range we want.
        chop_start_time = date_to_mirtime(ephem.Date(segments[i]['start_time'] - (cycle_time / 2) * ephem.second))
        chop_end_time = date_to_mirtime(ephem.Date(segments[i]['end_time'] + (cycle_time / 2) * ephem.second))
        chop_time_select = "time(%s,%s)" % (chop_start_time, chop_end_time)
        chop_file = "%s/segment_%04d.uvgen" % (args['--temp-dir'], i)
        if os.path.isdir(chop_file):
            shutil.rmtree(chop_file)
        #print "  Running uvaver to select time range %s - %s" % (chop_start_time, chop_end_time)
        with open("last_uvaver.dbg", "w") as fp:
            fp.write("uvaver vis=%s \"select=%s\" out=%s\n" % (simulated_name,
                                                               chop_time_select,
                                                               chop_file))
        miriad.uvaver(vis=simulated_name, select=chop_time_select,
                      out=chop_file)
        # Add this dataset to the uvcat string.
        uvcatString = add_uvcat(uvcatString, chop_file, args['--temp-dir'])

    # Do the final uvcat.
    finalout = "%s/allsegments.uvgen" % args['--temp-dir']
    print "  Concatenating final product"
    if os.path.isdir(finalout):
        shutil.rmtree(finalout)
    miriad.uvcat(vis=uvcatString, out=finalout)

    # Mix in the two datasets.
    print "  Adding the model to the initial dataset."
    if not os.path.isdir(finalout):
        print "Something has gone wrong and the fake dataset has not been generated."
        sys.exit()

    miriad.uvmodel(vis=args['<dataset>'], model=finalout, select="-auto", options="add",
                   out=args['--out'])
    print "Process complete. The mixed dataset can be found at %s" % args['--out']
    
if __name__ == '__main__':
    arguments = docopt(__doc__, version="Miriad Source Adding Helper 1.0")
    valid = True
    # Do some existence checking.
    if (('<dataset>' not in arguments) or
        (not os.path.isdir(arguments['<dataset>']))):
        print "No valid dataset found."
        valid = False
    if (('--out' not in arguments) and (valid == True)):
        # Make a default.
        arguments['--out'] = "%s.sourceadd" % arguments['<dataset>']
    # Check we have information to generate a source.
    if (('--source-file' not in arguments) or (arguments['--source-file'] is None)):
        if (('--ra' not in arguments) or ('--dec' not in arguments) or
            ('--flux' not in arguments) or ('--size' not in arguments) or
            ('--alpha' not in arguments)):
            print "Not enough information found to create a new source."
            valid = False
        elif ((len(arguments['--ra']) != len(arguments['--dec'])) or
              (len(arguments['--ra']) != len(arguments['--flux'])) or
              (len(arguments['--ra']) != len(arguments['--size'])) or
              (len(arguments['--ra']) != len(arguments['--alpha']))):
            print "Not every new source is properly specified."
            valid = False
    else:
        if (('--ra' in arguments) or ('--dec' in arguments) or
            ('--flux' in arguments) or ('--size' in arguments) or
            ('--alpha' in arguments)):
            print "A pre-made source file was specified, but so were new source parameters."
            valid = False
        elif (not os.path.isfile(arguments['--source-file'])):
            print "Specified source file cannot be found."
            valid = false
    # Check either the temp dir exists or can be made.
    if (not os.path.isdir(arguments['--temp-dir'])):
        os.makedirs(arguments['--temp-dir'])
    if (valid == True):
        print arguments
        add_source(arguments)

    
