<?php ini_set( "display_errors", 0); ?>
<html>
<head>
	<title>MOSAL - Multiobjective Sequence Alignment Tools</title>

	<!-- JQuery
	==================================================== -->	
	<script type='text/javascript' src="../js/jquery.min.js"></script>

	<!-- Bootstrap
	==================================================== -->
	<link href="../css/bootstrap.css" rel="stylesheet" type="text/css">
	<script type='text/javascript' src="../js/bootstrap-button.js"></script>
	<script type='text/javascript' src="../js/bootstrap.file-input.js"></script>
	<script type="text/javascript" src="../js/bootstrap-alert.js"></script>


	<meta name="description" content="Improvements on bicriteria pairwise sequence alignment: algorithms and applications" />
	<meta name="robots" content="index" />    
	<meta http-equiv="content-type" content="text/html;charset=UTF-8">

	<!-- Google Analytics ============================================== -->
	<script type="text/javascript">
		var _gaq = _gaq || [];
		_gaq.push(['_setAccount', 'UA-37516559-1']);
		_gaq.push(['_trackPageview']);

		(function() {
		var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
		ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
		var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
		})();
	</script>    
	<!-- End of Google Analytics ============================================== --> 

	<script type="text/javascript">
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

	<?php include '../includes/footer.html'; ?>

</body>
</html>
