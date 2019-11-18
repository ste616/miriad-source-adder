# miriad-source-adder
Helper script for adding a source into an existing Miriad data set

## Purpose
There are occasions when you'd like to add fake sources into a Miriad uv dataset, in order to test how
well they are found or characterised by some reduction method. Adding a point source to an existing dataset
is already possible using the task uvmodel, but it had been (until recently) impossible to add a more
complex source.

With the recent changes to Miriad, it is now possible to take one dataset and mix it with another, using uvmodel.
In that way, you could (in theory) model a complex source and create a dataset using uvgen, and mix that dataset with
the existing one with uvmodel.

This presents its own challenges however, as uvmodel needs the uv coordinates of the dataset containing the modelled
source to be the same as the existing dataset, since it only mixes them together; it cannot interpolate visibilities,
or generate what it should be given a complex model. The script presented in this repository handles the generation
of the new dataset using the exisiting dataset as a guide.

## Pre-requisites
This script requires that you have the latest ATNF Miriad installed, as changes to uvgen, uvmodel and the Miriad
library were made to support this functionality in November 2019.

Python modules required:
* docopt
* ephem
* mirpy
* numpy

## How it works
### Concept
For a single pointing dataset, this script is probably over-kill, although it will still work in the same way.
Imagine an observation with a 10 second cycle time, which pointed at a single spot in the sky, with coordinates
RA = 01:00:00, Dec = -20:00:00. Imagine then that you wanted to add a fake source at the coordinates
RA = 01:00:10, Dec = -20:59:50.

With uvgen as it currently stands, you would need to do the following:
* Determine the offset between the pointing location and the location where you want the fake source to go.
* Work out the starting and ending hour angles of this source on the date of the actual observation.
* Get the actual baseline lengths from the array which was used.
* Generate the fake dataset, and mix it with uvmodel.

None of this is particularly difficult to do, only time-consuming. Although, not everyone may be able to easily work
out how to get the actual baseline lengths. These actual lengths are crucial, as uvmodel can only mix datasets
together when they have exactly the same uv coordinates.

Things get quite difficult though when you're trying to mix in a fake source with a mosaicked dataset. In that case
you would need to do the same thing as above, for each individual pointing, and match the times of each pointing in the
generated dataset.

This script handles all that for you.

### In detail
Now let's look at how it does everything. For the moment though, we'll ignore how the script is actually run, but
instead just describe the steps the script makes.

__Step 1: Work out what was observed.__

This is easy: the script uses uvindex to make a list of each pointing centre, and when it was observed. The
script runs with a very short interval to ensure each cycle is listed.

__Step 2: Get the cycle time used during the observations.__

Again, easy, using uvlist this time with `options=variables,full`.

__Step 3: Get the positions of each antenna.__

Using uvlist with `options=array,full`.

__Step 4: Split the observing session into "segments".__

Each segment is a single source being observed for a contiguous block of time.

__Step 5: Generate a fake dataset for each segment.__

So, given the coordinates and parameters for your fake source, the script works out the inputs
needed for uvgen (offsets, hour angles, times etc.). It then runs uvgen and uvaver to generate
a dataset which precisely mirrors the real segment.

_With one caveat: during mosaicking, often the first cycle of a segment will be truncated at the start,
since the correlator will only begin recording data after the antennas get on source, which can be part-way
through a cycle. If one were to look at the uvindex outputs for the real vs. the fake dataset, this first
cycle may have a different start time, but the uv-coordinates are still taken from the mid-point of the
whole cycle, and thus are still perfectly compatible._

__Step 6: Concatenate all the segments together.__

The fake input to uvmodel needs to be the entire observation, where all segments are observed in the same
order and time as the real observation.

__Step 7: Mix the fake and real datasets together.__

Using uvmodel to do so.

## Using the script

Here's a fairly useful example command. We aim to add a fake Gaussian source at the coordinates:
RA = 22:45:01.6, Dec = -34:52:00. This source will have a major axis size of 10 arcsec, and a minor axis
size of 6 arcsec, at a position angle of 90 degrees. It will have a flux density of 800 mJy, and a spectral
index of -0.7.

Our real dataset is called real.5500, and we want to output a dataset called mixed.5500.

```
miriad-source-adder.py --ra 22:45:01.6 --dec -34:52:00 --size 10.0,6.0,90 -flux 0.8 --alpha -0.7 \
			--out mixed.5500 --temp-dir tmp real.5500
```

### Command line arguments

The script is controlled fully by the command line arguments you pass at runtime. The full list of arguments
and what they control is given below.

#### Default argument

You must always specify one argument, that being the name of the dataset to which you want to add a
source. This argument doesn't require any labelling. In the example above, `real.5500` is the default
argument.

#### Arguments without parameters

* `-h` or `--help`: show a brief usage guide listing the arguments supported and what they control.

#### Arguments requiring parameters

__Controlling the generated source__

To specify the parameters of the source to add to the dataset, use the following arguments.

* `-s` or `--source`, with parameter `FILE`: use the file `FILE` in the call to uvgen, passing it along as
the parameter to the argument `source`. While this may be completely appropriate if you're planning on adding
one or more sources to a single pointing dataset, it will not work (and this script will warn you) if you're planning
on adding a source to a mosaic. This is because the file expected by uvgen needs RA and Dec specified as offsets rather
than fixed positions on the sky.
* `-r` or `--ra`, with parameter `RA`: the right ascension of the source to add, given in sexagesimal (HH:MM:SS.S)
format. This should be in the same frame as the dataset to which you are adding the source (and will almost certainly
be J2000).
* `-d` or `--dec`, with parameter `DEC`: the declination of the source to add, given in sexagesimal (DDD:MM:SS.S)
format. Like RA, specify it in the same frame as the dataset.
* `-f` or `--flux`, with parameter `FLUX`: the flux density of the source to add, in Jy.
* `-S` or `--size`, with parameters `BMAJ,BMIN,BPA`: the size of the Gaussian source to add. The `BMAJ` and `BMIN`
parameters are the semi-major and semi-minor length of the source respectively, both in arcsec. The `BPA` is the position
angle of the Gaussian ellipsoid, in degrees, where 0 degrees puts the major axis along the right ascension direction,
and positive position angles rotate eastwise.
* `-a` or `--alpha`: the spectral index of the source to add.

__Controlling the code__

* `-o` or `--out`, with parameter `OUT`: the name of the output dataset, which is the combination of the input dataset
and the generated source.
* `-T` or `--temp-dir`, with parameter `DIR`: the name of a directory in which the script will put all the intermediate
datasets it needs to generate. After a successful execution, this directory can be safely deleted if desired, although the
script does not do so.


### Things you need to know

__You should apply the calibration on the real dataset beforehand__

The call to uvgen to create the fake dataset doesn't know anything about bandpass shape, or gains, or even system temperatures
that were present during the real observation. So we need to take the effect of all those things out beforehand.
Therefore, please do a full normal calibration of the real dataset, and apply it using uvaver, and pass the output of uvaver
as the input of this command.

__This script produces a LOT of files__

Depending on how many pointings there are in your mosaic, you will get a lot of new datasets from this script. There will be
one dataset for the initial generation, another for the clip to the proper time range, and then an output dataset each time
the script determines it needs to concatenate what it has together (there is a limit to the length of the string it can give
to uvcat for concatenation). There's also a text file to tell uvgen how to operate. For example, for the dataset I used to test
this script, which had 180 pointings, and 829 segments, a total of 1211 files were created.

This is why the script has a `--temp-dir` argument, and why you are encouraged to use it. All files will be generated under
this temporary directory. If you don't use it, all the files will be created in the current directory.

There is one exception: the final output dataset (as specified by the `--out` argument) is put in the current directory
(or in whichever directory you have specified as part of the `--out` argument).

__There is a default output name for the mixed dataset__

If you don't give an `--out` argument, then the mixed dataset will just be the name of the real dataset, with `.sourceadd`
appended.

__uvmodel doesn't always work__

Sometimes, uvmodel crashes with a memory error. In these cases you will need to run uvmodel yourself. The inputs to uvmodel are
(assuming the temporary directory is called tmp):

`uvmodel vis=real.5500 model=tmp/allsegments.uvgen select=-auto options=add out=mixed.5500`

The final concatenated dataset containing all the segments is always called `allsegments.uvgen`, under the temporary directory.

## Getting help

If you're having trouble with this script, please reach out to the owner of this repository (whose email address is not listed
here, for spam-mitigation reasons).
