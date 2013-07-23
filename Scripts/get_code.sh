#!/bin/sh
cp -r ../MultiObjective/ MOSAL
find MOSAL -name '.DS_Store' -exec rm -r {} \;
find MOSAL -name '.svn' -exec rm -r {} \;
find MOSAL -name 'EC' -exec rm -r {} \;
tar -zcvf MOSAL.tar.gz MOSAL
rm -r MOSAL