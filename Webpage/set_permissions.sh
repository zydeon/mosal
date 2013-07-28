#!/bin/sh
chmod -R 755 css/ img/ js/ plugins/ includes/
chmod 644 index.php
cd align;
chmod -R 755 MOSAL_code/
chmod -R 755 PHPMailer-master/
chmod 644 *.php
chmod 755 backend/
chmod 644 backend/.ht*
chmod 644 backend/*.php
chmod 666 backend/logs.html
chmod -R 777 inputs/
chmod 755 subscores
chmod 644 subscores/*
