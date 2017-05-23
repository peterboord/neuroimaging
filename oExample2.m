#!/usr/bin/octave -qf

arg_list = argv();
xFile = arg_list{1};

x = load("-ascii",xFile);

printf("%.4f", size(x));
