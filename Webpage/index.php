<?php ini_set( "display_errors", 0); ?>
<!DOCTYPE html>
<html lang="en">
<head>
	<?php $GLOBALS['path_to_root']=''; ?>
	<?php include 'includes/header.php'; ?>
</head>
<body>

	<div class="hero-unit">
		<h1>MOSAL <small>Multiobjective Sequence Alignment Tools</small></h1>
	</div>

	<div id="contents">
		<div class="row">
			<!--  Table of Contents   ================================================= -->
			<div class="span3" >
				<ul class="nav nav-tabs nav-stacked affix" style="left:40px;">
					<li><a href="#news">News</a></li>
					<li><a href="#description">Description</a></li>
					<li><a href="#license">License</a></li>
					<li><a href="#download">Download</a></li>
					<li><a href="#tryonline">Try on-line</a></li>
					<li><a href="#building">Building</a></li>
					<li><a href="#usage">Usage</a></li>
					<li><a href="#output">Output</a></li>
					<li><a href="#references">References</a></li>
					<li><a href="#acknowledgements">Acknowledgements</a></li>
				</ul>
			</div>
			<!--  End of Table of Contents   ================================================= -->


			<!--  Contents   ================================================= -->
			<div class="span9">
				<h2 id="news">News</h2>
				<ul>
					<li>Article "<a href="http://dx.doi.org/10.1186/1751-0473-9-2" target="_blank">MOSAL: Software tools for multiobjective sequence alignment</a>", accepted to <a href="http://www.scfbm.org" target="_blank">Source Code for Biology and Medicine</a>.</li>
					<li>BiMuSA, an implementation for biobjective multiple sequence alignment problems, is available <a href="http://www.zib.de/schenker/downloads/bimusa.tar.gz">here</a> (see <a href="http://eden.dei.uc.pt/~paquete/papers/MOSAL/R3.pdf
">technical report</a>)</li>
					<li>Due to an unexpected hardware failure, the web server was down several times in the last two weeks. The problem is now fixed.</li>
					<li>On-line version of MOSAL is available <a href="align">here</a></li>
					<li><a href="http://www.gcb2013.de" target="_blank">Poster accepted to the German Conference in Bionformatics 2013</a></li>
				</ul>
				<br>

				<h2 id="description">Description</h2>
				<p>MOSAL is a program for computing the Pareto optimal alignments for the bicriteria pairwise sequence alignment. 
			It allows the use of substitution matrices. </p>            
				<p>Support material for M.Abbasi, L. Paquete, A.Liefooghe, M. Pinheiro and P.Matias, 
				<a href="http://dx.doi.org/10.1093/bioinformatics/btt098">
			Improvements on bicriteria pairwise sequence alignment: algorithms and applications</a>, Bioinformatics, 29(8):996-1003, 2013.
		   </p>

				<p>Maintainer: <a href="mailto:mosal@dei.uc.pt">Pedro Matias <i class="icon-envelope"></i></a></p>

				<br>
				<h2 id="license">License</h2>

				  <p>This software is Copyright &copy; 2013 M.Abbasi, L. Paquete, A.Liefooghe, M. Pinheiro and P.Matias.</p>
				  <p>This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.</p>

				  <p>This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the <a href="http://www.gnu.org/licenses/gpl.html">GNU General Public License</a> for more details.</p>
				  <p>Appropriate reference to this software should be made when describing research in which it played a substantive role, so that it may be replicated and verified by others. The bicriteria sequence alignment problem and the algorithms which this software implements are described in detail in [<a href="#ref1">1</a>]. Please, mention the algorithm in the References section of your article. We would appreciate if you would email <a href="mailto:mosal@dei.uc.pt?Subject=Citation of paper referencing MOSAL" target="_top">mosal@dei.uc.pt</a> with citations of papers referencing this work.
				   MOSAL is distributed under the <a href="http://www.gnu.org/copyleft/gpl.html">GNU General Public License</a>. <p>

				<br>
				<h2 id="download">Download</h2>
				<p>Source code can be found <a href="MOSAL.tar.gz">here</a>.</p>

				<br>
				<h2 id="tryonline">Try on-line</h2>
				<p>A limited version of the application is available on-line <a href="align">here</a>.</p>            

				<br>
				<h2 id="building">Building</h2>
				<p>In GNU/Linux, the program can be compiled from source by invoking</p>
				
				<p><code>$ make</code></p>

				<br>
				<h2 id="usage">Usage</h2>

				<p><code>$ ./mosal seq1_file seq2_file [gaps|indels] [dp|dpp -b=NUMBER] [-ss=FILE] [--no-traceback]</code> </p>

				<table class="table">
					<thead>
						<tr>
							<th>Argument</th>
							<th>Description</th>
						</tr>            
					</thead>
					<tbody>
						<tr>
							<td><code>seq1_file</code></td>
							<td>path to the 1st sequence file (FASTA format)</td>
						</tr>
						<tr>
							<td><code>seq1_file</code></td>
							<td>path to the 2nd sequence file (FASTA format)</td>
						</tr>
						<tr>
							<td><code>gaps</code></td>
							<td>use gaps</td>
						</tr>
						<tr>
							<td><code>indels</code></td>
							<td>use indels</td>
						</tr>
						<tr>
							<td><code>dp</code></td>
							<td>use Dynamic Programming</td>
						</tr>  
						<tr>
							<td><code>dpp</code></td>
							<td>use Dynamic Programming with pruning (need to specify -b=NUMBER)</td>
						</tr>  
						<tr>
							<td><code>-ss=FILE</code></td>
							<td>use substitution score(ss) instead of #matches (path to the ss table)</td>
						</tr> 
						<tr>
							<td><code>-b=NUMBER</code></td>
							<td>specify the number of bounds for pruning version</td>
						</tr>
						<tr>
							<td><code>--no-traceback</code></td>
							<td>output only the scores without the alignments</td>
						</tr>   
					</tbody> 	    	   	    	    	
				</table>

				<br>
				<h2 id="output">Output</h2>
				<ul>
					<li>The non-traceback versions print only the score
					vectors. The first column corresponds to the number of
					matches/substitution scores and the second to the
					indels/gaps.</li>
					<li>The traceback versions print the scores and the
					  corresponding alignments</li>
				</ul>

				<br>
				<h2 id="references">References</h2>
				
				<table class="table">
					<tbody>
						<tr>
							<td class="span1" id="ref1">[1]</td>
							<td>M.Abbasi, L. Paquete, A.Liefooghe, M. Pinheiro and P.Matias, 
								<a href="http://dx.doi.org/10.1093/bioinformatics/btt098" target="_blank">Improvements on bicriteria pairwise sequence alignment: algorithms and applications</a>,
								Bioinformatics, 29(8):996-1003, 2013.
							</td>
						</tr>
						<tr>
							<td class="span1" id="ref2">[2]</td>
							<td>L. Paquete, P. Matias, M. Abbasi, M. Pinheiro
								<a href="http://dx.doi.org/10.1186/1751-0473-9-2" target="_blank">MOSAL: Software tools for multiobjective sequence alignment</a>,
								Source Code for Biology and Medicine, 9(2), 2014.
							</td>
						</tr>
					</tbody>
				</table>


				<br><br>
				<h2 id="acknowledgements">Acknowledgements</h2>
				This work was support by the Funda&ccedil;&atilde;o para a Ci&ecirc;ncia e
			Tecnologia, project MOSAL - Multiobjective sequence
				alignment (PTDC/EIA-CCO/098674/2008) and by FEDER, Programa Operacional
			Factores de Competitividade do QREN, ref. COMPETE: FCOMP-01-0124-FEDER-010024<p>
				
				<img width=15% src="img/Logo_Compete.gif"> 
				<img width=15% src="img/QREN_Logo(COR).png"> 
				<img width=15% src="img/UE_FEDER.jpg"> 
			
			</div>
		</div>
	</div>
	<!-- End of Contents   ================================================= -->

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

	<div class="navbar"  >
		<?php include 'includes/footer.php'; ?>
	</div>
</body>
</html>
