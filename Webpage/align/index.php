<?php ini_set( "display_errors", 0); ?>
<html>
<head>
    <?php $GLOBALS['path_to_root']='../'; ?>
	<?php include '../includes/header.php' ?>

	<script type="text/javascript">
    <?php $GLOBALS['path_to_root']=''; ?>
		$(document).ready(function() {
			// set events
			$("#alignm").click(function(){
				$('#frame').attr('src', 'input.php');
			});
			$("#visual").click(function(){
				$('#frame').attr('src', 'visual.html');
			});	
			$("#alignm").button('toggle');
		});	
	</script>

	<style type="text/css">
		iframe {
			width: 100%;
			height: 100%;
			border: none;
		}
	</style>	
</head>
<body>
	<br>
	<center>
		<div class="btn-group" data-toggle="buttons-radio">
			<input type="button" id="alignm" class="btn" value="Align sequences">
			<input type="button" id="visual" class="btn" value="Visualize results">
		</div>
		<hr>
	</center>

	<iframe id="frame" src="input.php"></iframe>

	<?php include '../includes/footer.php'; ?>

</body>
</html>
