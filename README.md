ProBez
======

Sometimes publications require 3D models to be shown.
Modern visualisation utilities mostly use raster graphics.
In the same time many models can be displayed as set of lines and Bezier curves.
The ProBez project intended to create software, which will approximate some 3D models (e.g. molecules) with lines and Bezier curves.

Requires
========

* CMake
* pkg-config
* GTK+-2.0 library and headers
* LibXml-2
* GSL
* [vmath](http://github.com/Lysovenko/vmath)
* [pgl](http://github.com/Lysovenko/pgl)

Installation
============

* Install dependencies if need.
* Make your build directory. E.g. <code>mkdir build</code> in folder with sources.
* Configure with CMake: <code>cmake ..</code> (replace <<..>> with path to sources.
* Do <code>make install</code> if you are using UNIX system or something like in other systems.

Using
=====

Wait couple of weeks while it become usable ;-)
The repo is supplied by an examples of a models.
Just open a model and drag mouse on the program's window.
