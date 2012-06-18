use strict;
use warnings;
use inc::Module::Install;

name 'PDL-FuncND';

license 'gpl3';
author 'Diab Jerius';

all_from 'lib/PDL/FuncND.pm';

requires 'PDL';

build_requires 'Test::More';

repository 'http://bitbucket.org/djerius/pdl-funcnd';
bugtracker 'https://rt.cpan.org/Public/Dist/Display.html?Name=PDL-FuncND';

author_tests('xt');

auto_install;

WriteAll;
