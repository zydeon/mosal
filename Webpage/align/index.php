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

	
	<title>Multiobjective Sequence Alignment</title>

	<!-- JQuery
	==================================================== -->	
	<script type='text/javascript' src="../js/jquery.min.js"></script>

	<!-- Bootstrap
	==================================================== -->
	<link href="../css/bootstrap.css" rel="stylesheet" type="text/css">
	<script type='text/javascript' src="../js/bootstrap-button.js"></script>
	<script type='text/javascript' src="../js/bootstrap.file-input.js"></script>
	<script type="text/javascript" src="../js/bootstrap-alert.js"></script>

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
	<script type="text/javascript" src="../js/align.js"></script>
	<script type="text/javascript" src="../js/visualize.js"></script>
	<link href="../css/mystyles.css" rel="stylesheet" type="text/css">

</head>
<body>
	<center>
		<h2>Multiobjective Sequence Alignment</h2>
		<br>
		<div id="alerts_container"></div>
		<form action="javascript:startAlignment()">
			<fieldset class="step">
				<legend class="step">Step 1 - Sequences input<br><small>(size &lt; 2000)</small></legend>
				<div id="alerts_seq1"></div>
				<strong>Sequence 1</strong><br>
				<textarea id="seq1" maxlength="2000" placeholder="ATGAACAATCAAGCATACGGTGTTACACC..." style="width:75%;" rows="4" required></textarea><br>
				<br>
				<div id="alerts_seq2"></div>
				<strong>Sequence 2</strong><br>
				<textarea id="seq2" maxlength="2000" placeholder="CTGACCACGAAGACATACGGAGTAACTGA..." style="width:75%;" rows="4" required></textarea><br>
			</fieldset>

			<fieldset class="step">
				<legend class="step">Step 2 - Score function</legend>
				<div class="btn-group" data-toggle="buttons-radxio">
					<input type="button" id="ibtn" class="btn" value="Indels">
					<input type="button" id="gbtn" class="btn" value="Gaps">
				</div>
				<br>
				<br>
				<div class="btn-group" data-toggle="buttons-radio">
					<input type="button" id="mbtn" class="btn" value="Matches">
					<input type="button" id="sbtn" class="btn" value="Substitution score">
				</div>
				<br>
				<div id="subscore_form">
					<br>
					<div id="alerts_subscore"></div>
					<strong>Choose a substitution score matrix</strong><br>
					<!-- fetch substitution score matrixes -->
					<?php
						$dir = opendir("subscores");
						while($entry = readdir($dir))
							if(preg_match("/^BLOSUM.*\.tbl$/", $entry))
								$blosums[] = substr($entry, 0, strrpos($entry, '.'));
							else if(preg_match("/^PAM.*\.tbl$/", $entry))
								$pams[] = substr($entry, 0, strrpos($entry, '.'));
						closedir($dir);
					?>
					<select id="subscore_select">
						<option value="-">-</option>
						<optgroup label="BLOSUM">
							<?php foreach ($blosums as $b) echo "<option value='$b'>$b</option>"; ?>
						</optgroup>
						<optgroup label="PAM">
							<?php foreach ($pams as $p) echo "<option value='$p'>$p</option>"; ?>
						</optgroup>
					</select>
					<p>or</p>
					<label for="subscore_upl"><i class="icon-circle-arrow-up"></i><strong> Upload substitution score matrix</strong></label>
					<input type="file" id="subscore_upl">
				</div>
			</fieldset>

			<fieldset class="step">
				<legend class="step">Step 3 - Sequence alignment options</legend>
				<table>
					<tr>
						<td><input type="checkbox" id="tb_chkbox"></td>
						<td>Traceback (display alignments)</td>
					</tr>					
					<tr>
						<td><input type="checkbox" id="pruning_chkbox"></td>
						<td>Use pruning</td>
					</tr>
				</table>
				<div id="pruning_form">
					<br>							
					<label for="pruning"><strong>Choose number of bounds</strong></label>
					<input type="range" min="1" max="25" id="pruning">
					<input type="text" id="pruningValue" maxlength="2">							
				</div>

			</fieldset>

			<fieldset class="step">
				<legend class="step">Step 4 - Submit</legend>
				<input type="checkbox" id="email_chkbox">Send results to e-mail
				<br>
				<br>
				<div id="email_form"></div>
				<input type="submit" class="btn btn-primary" value="Align">
			</fieldset>
		</form>

		<hr>
		<div id="load"></div>
		<div id="results"></div>
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

	</center>
</body>
</html>
