<?php

function add_to_logs($msg, $cmd="[no command]"){
	$host = $_SERVER['REMOTE_ADDR'];
	$port = $_SERVER['REMOTE_PORT'];
	date_default_timezone_set('Europe/Lisbon');
	$date = date('m/d/Y h:i:s a');

	$f = fopen('backend/logs.html', 'a');
	if($f){
		fwrite($f,
					"<table rules='rows'>".
						"<tr>".
							"<td class='head'>DATE</td>".
							"<td>$date</td>".
						"</tr>".
						"<tr>".
							"<td class='head'>USER HOST</td>".
							"<td>$host</td>".
						"</tr>".
						"<tr>".
							"<td class='head'>USER PORT</td>".
							"<td>$port</td>".
						"</tr>".
						"<tr>".
							"<td class='head'>MSG</td>".
							"<td>$msg</td>".
						"</tr>".
						"<tr>".
							"<td class='head'>CMD EXECUTED</td>".
							"<td>$cmd</td>".
						"</tr>".
						"<br>".
						"<br>".
					"</table>"
		);

		fclose($f);
	}
}

function log_($text){
	echo "<script>console.log('$text')</script>";
}

?>

