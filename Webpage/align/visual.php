<?php ini_set( "display_errors", 0); ?>
<html>
<head>
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
	<title>Visualize multiobjective sequence alignments</title>

	<!-- JQuery
	==================================================== -->	
	<script type='text/javascript' src="../js/jquery.min.js"></script>

	<!-- Bootstrap
	==================================================== -->
	<link href="../css/bootstrap.css" rel="stylesheet" type="text/css">
	<script type='text/javascript' src="../js/bootstrap-button.js"></script>
	<script type='text/javascript' src="../js/bootstrap.file-input.js"></script>

	<!-- JQPlot
	==================================================== -->	
	<link href="../css/jquery.jqplot.min.css" rel="stylesheet" type="text/css">
	<script type='text/javascript' src="../js/jquery.jqplot.min.js"></script>
	<script type="text/javascript" src="../plugins/jqplot.highlighter.min.js"></script>
	<script type="text/javascript" src="../plugins/jqplot.canvasTextRenderer.min.js"></script>
	<script type="text/javascript" src="../plugins/jqplot.canvasAxisLabelRenderer.min.js"></script>
	<script type="text/javascript" src="../plugins/jqplot.cursor.min.js"></script>

	<!-- My includes
	==================================================== -->	
	<script type="text/javascript" src="../js/visualize.js"></script>
	<script type="text/javascript" src="../js/upload.js"></script>
	<script type="text/javascript" src="../js/bootstrap-alert.js"></script>
	<link href="../css/mystyles.css" rel="stylesheet" type="text/css">
	<script type="text/javascript">
		$(document).ready(function() {
			$('input[type=file]').bootstrapFileInput();
		});
	</script>

</head>
<body>
	<center>
		<h2>Visualize multiobjective sequence alignments</h2>
		<br>
		<div id="alerts_container"></div>
		<form action="javascript:executeUpload()">
			<label for="data_file"><i class="icon-circle-arrow-up"></i><strong> Upload either values or alignments file</strong></label>
			<input type="file" accept="text/tab-separated-values,text/plain,text/csv,text/html,text/richtext,text/rtf,text/xml" id="data_file" required>
			<br>
			<br>
			<input type="submit" class="btn btn-primary" value="Visualize">
		</form>
			<hr>
			<div id="results"></div>

			<!-- Leave space to footer -->
			<br>
			<br>
			<br>
			<br>
			<br>
			<br>
			<br>
			<br>
			<br>
			<br>
			<br>
			<!-- End of Leave space to footer -->
	
			<div class="navbar">
				<?php include '../includes/footer.php'; ?>
			</div>
	</center>
</body>
</html>
