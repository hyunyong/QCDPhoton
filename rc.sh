#!/bin/bash
g++ calN.cpp -o calN.o `root-config --cflags --glibs`
./calN.o
