<?php ini_set( "display_errors", 0); ?>
<html>
<head>
<?php $GLOBALS['path_to_root']='../'; ?>
<?php include '../includes/header.php' ?>

	<script type="text/javascript">
		$(document).ready(function() {
			// set events
			$("#alignm").click(function(){
				$('#frame').attr('src', 'input.php');
			});
			$("#visual").click(function(){
				$('#frame').attr('src', 'visual.php');
			});
			$("#help").click(function(){
				$('#frame').attr('src', 'help.php');
			});
			$("#alignm").button('toggle');
		});
	</script>

	<style type="text/css">
		body{
			overflow: hidden;
		}
		iframe {
			width: 100%;
			height: 90%;
			border: none;
		}
		#header{
			/*background-color: #EEEEEE;*/
			padding-top: 1em;
		}
	</style>
</head>
<body>
	<div id="header">
		<center>
			<div class="btn-group" data-toggle="buttons-radio">
				<input type="button" id="alignm" class="btn" value="Align sequences">
				<input type="button" id="visual" class="btn" value="Visualize results">
				<input type="button" id="help"   class="btn" value="Help">
			</div>
		<hr>
		</center>
	</div>

	<iframe id="frame" src="input.php"></iframe>
</body>
</html>
