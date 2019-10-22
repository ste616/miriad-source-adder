#!/usr/bin/python
"""Miriad Source Adding Helper

Usage:
  miriad-source-adder.py [--source-file=<str>] [--ra=<str> ...] [--dec=<str> ...] [--flux=<flux> ...] [--size=<bmaj,bmin,bpa> ...] [--alpha=<alpha> ...] [--out=<str>] <dataset>

-h --help               show this
-s --source-file FILE   file to give to uvgen as "source" parameter
-r --ra RA              the RA of a source to add HH:MM:SS.S
-d --dec DEC            the declination of a source to add DDD:MM:SS.S
-f --flux FLUX          the flux density of a source to add in Jy
-S --size BMAJ,BMIN,BPA the size of a source to add, in arcsec for each number
-a --alpha ALPHA        the spectral index of a source to add
-o --out OUT            output the mixed dataset with this name

"""

from docopt import docopt
import math
from mirpy import miriad
import ephem
from datetime import date, datetime
import os
import shutil
import numpy as np

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

def add_source(args):
    # Let's create a uvindex of the dataset first.
    dataset_name = args['<dataset>']
    miriad.set_filter('uvindex', filter_uvindex)
    index_data = miriad.uvindex(vis=dataset_name, interval="0.1")

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
    
    for i in xrange(0, len(index_data['sources'])):
        # Make a source object for this position.
        source = ephem.FixedBody()
        source.name = index_data['sources'][i]['name']
        source._epoch = "2000"
        source._ra = index_data['sources'][i]['right_ascension']
        source._dec = index_data['sources'][i]['declination']

        # Find the times that this source was observed.
        source_observation_idx = np.where(index_data['index']['source'] == source.name)
        source_observation_times = index_data['index']['time'][source_observation_idx]
        
        # Where was the source at the start of the observations?
        # Go to the starting time of the observation.
        telescope.date = ephem.Date(source_observation_times[0])
        # We really just need the hour angle.
        start_hour_angle = get_hour_angle(source, telescope)

        # And the same for the end of the observation.
        telescope.date = ephem.Date(source_observation_times[-1])
        finish_hour_angle = get_hour_angle(source, telescope)
        print "%s %f %f" % (source.name, start_hour_angle, finish_hour_angle)

        # Work out the 0 hour angle. Get the nearest time to 0.
        transit_time = None
        if (abs(start_hour_angle) < abs(finish_hour_angle)):
            transit_time = source_observation_times[0] - (start_hour_angle * ephem.hour)
        else:
            transit_time = source_observation_times[-1] - (finish_hour_angle * ephem.hour)
        transit_time = ephem.Date(transit_time)

        print "transit time is %s" % transit_time
        rounded_transit_time = round_date_5second(transit_time)
        print "rounded transit time is %s" % rounded_transit_time
        telescope.date = transit_time
        print "LST at transit time (check) is %s" % telescope.sidereal_time()

        # Work out the inputs to uvgen.
        if (arguments['--source-file'] is not None):
            uvgen_source = arguments['--source-file']
        # Determine the offset between the sources we want to add and the pointing centre.
        for j in xrange(0, len(args['ra'])):
            nsource = ephem.FixedBody()
            nsource._epoch = "2000"
            nsource._ra = args['ra'][j]
            nsource._dec = args['dec'][j]
        
        uvgen_source = "addsource" # FIX THIS HARDCODED VALUE
        uvgen_telescop = telescope_coordinates['telescope'].lower()
        uvgen_stokes = ",".join(index_data['polarisations']).lower()
        print uvgen_source
        print uvgen_telescop
        print uvgen_stokes
        uvgen_lat = "%.6f" % math.degrees(telescope.lat)
        print uvgen_lat
        uvgen_radec = "%s,%s" % (index_data['sources'][0]['right_ascension'],
                                 index_data['sources'][0]['declination'])
        print uvgen_radec
        # Extend the HA range a little on each side
        uvgen_harange = "%d,%d,%.16f" % (math.floor(start_hour_angle),
                                         math.ceil(finish_hour_angle), (10. * 1.00273790935 / 3600.))
        print uvgen_harange
        # The time at 0 HA.
        uvgen_time = date_to_mirtime(rounded_transit_time)
        print uvgen_time
        # Use only frequency config 0, but we could probably change that if we need to.
        print index_data['freq_configs']
        chan_offset = math.floor(index_data['freq_configs'][0]['nchannels'][0] / 2)
        chan_spacing_mhz = index_data['freq_configs'][0]['frequency_increment'][0] * 1000.
        freq_offset_mhz = chan_spacing_mhz * chan_offset
        width_mhz = chan_spacing_mhz * index_data['freq_configs'][0]['nchannels'][0]
        uvgen_corr = "%d,1,%d,%.3f" % (index_data['freq_configs'][0]['nchannels'][0],
                                       freq_offset_mhz, width_mhz)
        print uvgen_corr
        uvgen_freq = "%.3f,0.0" % index_data['freq_configs'][0]['frequency1'][0]
        print uvgen_freq
        # Now we have to write out the antenna locations file.
        uvgen_ant = "antenna_configuration.file"
        output_antenna_file(telescope_coordinates, uvgen_ant)
        uvgen_baseunit = 3.33564
        
        # The name of the simulated dataset.
        simulated_name = "%s.uvgen" % dataset_name
        # Delete this set if it exists.
        if os.path.isdir(simulated_name):
            shutil.rmtree(simulated_name)

        # Run uvgen now.
        miriad.uvgen(source=uvgen_source, ant=uvgen_ant, baseunit=uvgen_baseunit,
                     telescop=uvgen_telescop, corr=uvgen_corr, time=uvgen_time,
                     freq=uvgen_freq, radec=uvgen_radec, harange=uvgen_harange,
                     stokes=uvgen_stokes, lat=uvgen_lat, out=simulated_name)


if __name__ == '__main__':
    arguments = docopt(__doc__, version="SourceAdd 1.0")
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
    if (valid == True):
        addsource(arguments)

    
