#!/bin/bash

./externals_all << endofinput > runout/$1.out
$1
endofinput
