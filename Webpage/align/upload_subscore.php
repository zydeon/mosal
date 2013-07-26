<?php
	include 'utils.php';

	if ($_FILES["subscore_file"]){
		if ($_FILES["subscore_file"]["error"] > 0){
			$response['status'] = 'error';
			add_to_logs("Error on uploading substitution score file.");
		}
		else {
			$timestamp = $_POST["timestamp"];
			$input_path = 'inputs/'.$timestamp.'-'.$_SERVER['REMOTE_ADDR'].'-'.$_SERVER['REMOTE_PORT'];

			if(!file_exists($input_path))
				mkdir($input_path);

			// instead of move_uploaded_file (avoid dangerous file extensions)
			$file_content = file_get_contents($_FILES["subscore_file"]["tmp_name"]);
			$f = fopen($input_path.'/subscore.tbl', 'w');
			fwrite($f, $file_content);
			fclose($f);

			$response['status'] = 'success';
		}		
	}

	echo json_encode($response);
?>
