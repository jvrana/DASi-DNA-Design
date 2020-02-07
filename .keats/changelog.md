# dasi change log
## 0.0.19a
**2020-02-07T07:09:29.642357**





## 0.0.19
**2020-02-06T18:28:08.818234**
bug fixes

 - bug fixes for scoring misprimings


## 0.0.18
**2020-02-06T13:52:04.675546**
bug fix

 - fixes serious bug where sequence scoring via `DNAStats` was case sensitive to the sequence


## 0.0.17
**2020-02-06T13:08:56.007142**
bug fixes

 - fixes bugs in scoring misprimings


## 0.0.16
**2020-02-06T12:48:09.054983**
speed improvements

 - scoring misprimings is much faster
 - exposed post processing parameters to `design.compile`


## 0.0.15
**2020-02-05T17:22:51.226522**
features

 - scores PCR products


## 0.0.14
**2020-02-03T08:21:23.810686**
bug fixes

 - many bug fixes to pass tests


## 0.0.13
**2020-02-03T07:31:02.036588**
bug fixes




## 0.0.12
**2020-02-02T23:21:51.948271**
new featues; output JSON format

 - New validated JSON formate for design outputs. Accessible by `design.out()`
 - JSON Schemas are grouped in `dasi.schemas`
 - updated documentation
 - improved testing
 - minor bug fixes
 - `design.report()` pulls up a new report. Plot coverage using `design.report().plot_coverage(show=True)`
 - `optimize_library` and `compile_library` have been removed from LibraryDesign. Use `optimize` and `compile` instead.


## 0.0.11
**2020-01-29T07:42:52.819562**
feature

 - cached span cost in dasi.cost module


## 0.0.10
**2020-01-28T14:02:29.086775**
bug fixes for failed builds




## 0.0.9
**2020-01-09T15:58:22.196391**
bug fixes for gaps on edge of region




## 0.0.8
**2020-01-09T14:32:35.169279**
library design support




## 0.0.7a0
**2019-11-27T09:44:03.094130**
features for LibraryDesign




## 0.0.6
**2019-10-18T12:49:46.896456**
update API breaking dependencies

 - update networkx to 2.4, update pyblast to 0.5


## 0.0.5
**2019-10-17T12:14:53.728652**
setup updates




## 0.0.4
**2019-10-17T11:47:06.003351**
command-line interface




## 0.0.2
**2019-10-01T15:25:54.613524**





## 0.1.0
**2019-07-13T10:01:24.258900**
project initialized


