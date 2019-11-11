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

## How it works
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
