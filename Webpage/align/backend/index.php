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
</head>
<body>
	<h1>Logs</h1>
	<?php echo file_get_contents('logs.html'); ?>
</body>
</html>