<?php
	ini_set( "display_errors", 0);

	$f = fopen('logs.html', 'w');
	if($f){
		fwrite($f, '');
		fclose($f);

		$response['status'] = 'success';
	}
	else{
		$response['status'] = 'error';
		$response['debug'] = 'Could not open logs file!';
	}

	echo json_encode($response);
?>