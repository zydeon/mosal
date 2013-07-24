<?php

function add_to_logs($msg, $cmd="[no command]"){
	$host = $_SERVER['REMOTE_ADDR'];
	$port = $_SERVER['REMOTE_PORT'];
	date_default_timezone_set('Europe/Lisbon');
	$date = date('m/d/Y h:i:s a');

	$f = fopen('backend/logs.txt', 'a');
	if($f){
		fwrite($f, "DATE\t$date\nHOST\t$host\nPORT\t$port\nMSG\t$msg\nCMD\t$cmd\n########################\n");
		fclose($f);
	}
}


?>