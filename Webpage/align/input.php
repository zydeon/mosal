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
	<style type="text/css">
		.nav{
			/*width:200px;*/
			background-color:white;
			left:20px;
			top:100px;
		}
		textarea{
			width:75%;
			margin-bottom:0px;
			font-family: Monaco, Menlo, Consolas, "Courier New", monospace;
		}
	</style>

</head>
<body>
	<center>
		<h2>Multiobjective Sequence Alignment</h2>
	</center>
		<br>

			<ul class="nav nav-tabs nav-stacked affix">
				<li><a href="#step1">Step 1 - Sequences input</a></li>
				<li><a href="#step2">Step 2 - Score function</a></li>
				<li><a href="#step3">Step 3 - Alignment options</a></li>
				<li><a href="#step4">Step 4 - Submit</a></li>
			</ul>

		<center>
		<!--  Contents   ================================================= -->
		<div id="alerts_container"></div>
		<form action="javascript:startAlignment()">
			<fieldset id="step1" class="step">
				<legend class="step">Step 1 - Sequences input (size &lt; 2000)<br><small><strong>FASTA format</strong></small></legend>
				<div class="btn-group" data-toggle="buttons-radio">
					<input type="button" id="protein" class="btn" value="Protein">
					<input type="button" id="dna" class="btn" value="DNA">
				</div>
				<br>
				<br>
				<div id="alerts_seq1"></div>
				<strong>Sequence 1</strong><br>
				<textarea id="seq1" maxlength="2000" rows="4" required></textarea>
				<br>
				<button id="sample1" type="button" class="btn btn-link">sample</button>
				<br>
				<br> 
				<div id="alerts_seq2"></div>
				<strong>Sequence 2</strong><br>
				<textarea id="seq2" maxlength="2000" rows="4" required></textarea>
				<br>
				<button id="sample2" type="button" class="btn btn-link">sample</button>				
			</fieldset>

			<fieldset id="step2" class="step">
				<legend class="step">Step 2 - Score function</legend>
				<div class="btn-group" data-toggle="buttons-radio">
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
					<?php displaySelect("Protein", "protein_select") ?>
					<?php displaySelect("DNA", "dna_select", true) ?>

					<p>or</p>
					<label for="subscore_upl"><i class="icon-circle-arrow-up"></i><strong> Upload substitution score matrix</strong></label>
					<input type="file" id="subscore_upl">
				</div>
			</fieldset>

			<fieldset id="step3" class="step">
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

			<fieldset id="step4" class="step">
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
		<div class="navbar">
			<?php include '../includes/footer.php'; ?>
		</div>
	</center>
</body>
</html>

<?php
	function getSubscoreTables($path){
		$groups['undefined'] = array();
		$dir = opendir($path);
		while($entry = readdir($dir)){
			if(preg_match("/^[^\.]+\.?[^\.]+$/", $entry)){		// hide hidden files
				if(is_dir($path."/".$entry)) 	$groups[$entry] = array();
				else 							$groups['undefined'][] = $entry;
			}
		}
		closedir($dir);

		/* read groups (sub-folders) */
		foreach ($groups as $g => $tables) {
			if($g != 'undefined'){
				$dir = opendir($path."/".$g);
				while($entry = readdir($dir)){
					if(preg_match("/^[^\.]+\.?[^\.]+$/", $entry)){			// hide hidden files
						$groups[$g][] = $entry;
					}
				}
				closedir($dir);
			}
		}
		return $groups;
	}

	function displaySelect($type, $id, $hidden=false){
		$groups = getSubscoreTables("subscores/$type");
		$display = $hidden==true ? "style='display:none;'":"";
		echo "<select id='$id' $display><option value='-'>-</option>";
		
		foreach ($groups['undefined'] as $t) {
			echo "<option value='$type/$t'>$t</option>";
		}
		foreach ($groups as $g => $tables) {
			if($g!='undefined'){
				echo "<optgroup label='$g'>";
				foreach ($tables as $t) echo "<option value='$type/$g/$t'>".stripExt($t)."</option>";
				echo "</optgroup>";
			}
		}
		echo "</select>";		
	}

	function stripExt($filename){
		return substr($filename, 0, strrpos($filename, '.'));
	}
?>	
