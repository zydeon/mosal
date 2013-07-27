<?php ini_set( "display_errors", 0); ?>


<html>
<head>
	<title>Backend</title>
	<!-- JQuery
	==================================================== -->	
	<script type='text/javascript' src="../../js/jquery.min.js"></script>
	<style type="text/css">
		iframe {
			width: 100%;
			height: 100%;
			border: none;
		}
		td.head{
			text-align:right;
			font-weight: bold;
			padding: 0em 1em 0em 0em;
		}
	</style>
	<script type="text/javascript">
		function deleteLogs(){
			var sure = confirm("Do you really want to delete all logs ?");
			if(sure ==true){
				jQuery.get('delete_logs.php', function(response){
					var response = JSON.parse(response);
					if(response.status=='success')
						window.location.href = window.location.href; 		// reload page
					else
						alert(response.debug);
				})
			}
		}
	</script>
</head>
<body>
	<h1>Logs</h1>
	<button onclick="deleteLogs();">Delete all logs</button>
	<?php echo file_get_contents('logs.html'); ?>
</body>
</html>