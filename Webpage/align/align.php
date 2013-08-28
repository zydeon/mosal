<?php
	ini_set( "display_errors", 0);

	include 'utils.php';
	require 'PHPMailer-master/class.phpmailer.php';

	$timestamp = $_POST['timestamp'];
	$input_path = 'inputs/'.$timestamp.'-'.$_SERVER['REMOTE_ADDR'].'-'.$_SERVER['REMOTE_PORT'];
	$executable_path = '../../MultiObjective';
	$remove_temporary_folder = true;

	$output    = array();
	$seq1      = $_POST['seq1'];
	$seq2      = $_POST['seq2'];
	$problem   = $_POST['is_indels']	== "true" 	? "indels" : "gaps";
	$approach  = $_POST['bound_number']	== "0" 		? "dp" 	: "dpp -b=".$_POST['bound_number']; 
	$traceback = $_POST['is_traceback']	== "true" 	? "" 		: "--no-traceback";
	$email     = $_POST['email'];
	$subject   = $_POST['email_subject'];
	$subscore  = setSubscore($_POST['subscore'], $input_path);

	// create temporary files
	if(!file_exists($input_path))
		mkdir($input_path);	
	// sequence 1
	$f = fopen($input_path.'/s1.fasta', 'w');
	fwrite($f, format_seq($seq1));
	fclose($f);	
	// sequence 2
	$f = fopen($input_path.'/s2.fasta', 'w');
	fwrite($f, format_seq($seq2));
	fclose($f);		

	$cmd = $executable_path."/mosal ".$input_path."/s1.fasta ".$input_path."/s2.fasta ".$problem." ".$approach." ".$subscore." ".$traceback;
	exec($cmd, $output, $return );
	if($return == 0){
		// success
		$response['status'] = 'success';
		$response['data']['values'] = array();

		if($_POST['is_traceback'] == "true"){
			$response['data']['alignments'] = array();

			for ($i=0; $i < sizeof($output); $i+=4) {
				// get solution values
				sscanf($output[$i], "%*[^=]=%d %*[^=]=%d", $a, $b);
				$response['data']['values'][] = array($a, $b);

				// get alignments
				$response['data']['alignments'][] = array($output[$i+1], $output[$i+2]);
			}
		}
		else{
			foreach ($output as $score) {							// format as array of arrays [[a1 b2], [a2 b2],...]
				sscanf($score, "%d %d", $a, $b);
				$response['data']['values'][] = array($a, $b);
			}			
		}
		// $response['debug'] = json_encode($cmd);

		// send email
		if($email != ""){
			$mail = new PHPMailer;

			$mail->IsSMTP();
			$mail->Host = 'smtp.dei.uc.pt';  

			$mail->From = 'mosal@dei.uc.pt';
			$mail->FromName = 'MOSAL webmaster';
			$mail->AddAddress($email);

			$mail->AddAttachment($input_path.'/s1.fasta');
			$mail->AddAttachment($input_path.'/s2.fasta');

			// create output files
			$s1 = $_POST['subscore']	== "true" 	? "subscore" : "matches";
			$s2 = $_POST['is_indels']	== "true" 	? "indels" : "gaps";
			$header = $s1 . "\t" . $s2;

			$fp = fopen($input_path.'/values.tsv', 'w');
			fwrite($fp, output_values($header, $response['data']['values']));
			fclose($fp);
			$mail->AddAttachment($input_path.'/values.tsv');
			if($_POST['is_traceback'] == "true"){
				$fp = fopen($input_path.'/alignments.txt', 'w');
				fwrite($fp, output_aligns($header, $response['data']['values'], $response['data']['alignments']) );
				fclose($fp);
				$mail->AddAttachment($input_path.'/alignments.txt');
			}

			$mail->Subject = $subject;
			$mail->Body = 'In the attachments section you will find the results of the alignment solutions along with the sequences you requested.'."\n";
			$mail->Body .= 'Go to http://mosal.dei.uc.pt/align/ so you can upload the results and visualize them.';

			if(!$mail->Send()){
				$msg = "Could not send email to $email!\n";
				$msg .= 'Mailer Error: ' . $mail->ErrorInfo;
				$response['debug'] = json_encode($msg);
				// $remove_temporary_folder = false;
				add_to_logs($msg, $cmd);
			}
		}

		if($remove_temporary_folder){
			// remove temporary folder
			exec("rm -r $input_path", $output, $return);
			if($return != 0){
				$response['debug'] = json_encode($output);
				add_to_logs("Error on removing inputs temporary folder.", $cmd);
			}
		}
	}
	else{
		$response['status'] = 'error';
		$response['debug'] = "Error on executing the program";
		$response['debug'] = json_encode($cmd);
		add_to_logs("Error on running the executable.", $cmd);
	}
	echo json_encode($response);


	function format_seq($s){
		if(!isFASTA($s)){
			// create fasta files content from string
			$inc = 70;
			$s2 = "> temporary file of sequence input\r\n";
			for ($i=0; $i < strlen($s); $i+=$inc ) 
				$s2.= substr($s, $i, min($inc, strlen($s)-$i))."\r\n";

			return $s2;
		}
		return $s;
	}
	function output_values($header, $values){
		// return output of data values to be stored in a file
		$contents = $header . "\r\n";
		foreach ($values as $key => $value)
			$contents .= $value[0] . "\t" . $value[1] . "\r\n";
		return $contents;
	}
	function output_aligns($header, $values, $alignments){
		$contents = $header . "\r\n";
		foreach ($values as $index => $value) 
			$contents .= $value[0] . "\t" . $value[1] . "\r\n" .
						 $alignments[$index][0] . "\r\n" .
						 $alignments[$index][1] . "\r\n";
		
		return $contents;
	}
	function setSubscore($formdata, $input_path){
		switch($formdata){
			case "false": 	return "";
			case "upload": 	return "-ss=\"".$input_path."/subscore.tbl\"";
			default: 		return "-ss=\"subscores/".$formdata."\"";
		}
	}
	function isFASTA($seq){
		return preg_match("/^>.*(\r?\n[A-Z\*]+)+$/", $seq);
	}
?>
