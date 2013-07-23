<?php

	if ($_FILES["subscore_file"]){
		if ($_FILES["subscore_file"]["error"] > 0){
			$response['status'] = 'error';
		}
		else {
			$timestamp = $_POST["timestamp"];

			if(!file_exists("inputs/".$timestamp))
				mkdir("inputs/".$timestamp);

			// instead of move_uploaded_file (avoid dangerous file extensions)
			$file_content = file_get_contents($_FILES["subscore_file"]["tmp_name"]);
			$f = fopen('inputs/'.$timestamp.'/subscore.tbl', 'w');
			fwrite($f, $file_content);
			fclose($f);

			$response['status'] = 'success';
		}		
	}

	echo json_encode($response);
?>
