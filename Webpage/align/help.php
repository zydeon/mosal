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

    <!-- My includes
    ==================================================== -->    
    <link href="../css/mystyles.css" rel="stylesheet" type="text/css">
    <style type="text/css">
        body{
            margin-left:25%;
            margin-right:25%;
        }
        .nav{
            /*width:200px;*/
            background-color:white;
            left:15%;
            top:8%;
        }
    </style>

</head>
<body>

        <ul class="nav nav-tabs nav-stacked affix">
            <li><a href="#step1">Step 1</a></li>
            <li><a href="#step2">Step 2</a></li>
            <li><a href="#step3">Step 3</a></li>
            <li><a href="#step4">Step 4</a></li>
        </ul>    

        <br>
        <h2>Help</h2> <br>
        <p>This webpage describes the process of using <a href="http://mosal.dei.uc.pt/align" target="_parent">MOSAL</a> to get multiobjective sequence alignments. In order to use the tool, one must follow a multiple step process that is explained below.</p>
        <br>
        <h3 id="step1">Step 1 - Sequences input</h3> <br>
        <p>This is the step where the user provide the target sequences. The format of the sequences can either be the <a href="http://en.wikipedia.org/wiki/FASTA_format" target="_blank">FASTA</a> format or just a string of characters. In both formats only upper-case letters are allowed, as well as the "*" symbol. </p>
        <p>The selection between <b>Protein</b> and <b>DNA</b> influences the choice of the substitution score matrix (step 2). The buttons below each sequence input field allow to introduce an example of a FASTA formatted sequence.</p>

        <br>
        <h3 id="step2">Step 2 - Score function</h3> <br>
        <p>In this step the user must choose which objectives to consider:</p>
        <ul>
            <li>Minimize either <b>Indels</b> or <b>Gaps</b></li>
            <li>Maximize either the number of <b>Matches</b> or the <b>Substitution score</b> </li>
        </ul>
        <p>The user can choose one of the default substitution score matrixes (dependent on choice between Protein or DNA) or upload one of its own as long as it follows the standards (see <a href="http://en.wikipedia.org/wiki/Substitution_matrix" target="_blank">link</a>).</p>

        <br>
        <h3 id="step3">Step 3 - Alignment options</h3> <br>
        <p>The "traceback" option in this step, when selected, allow to compute both the sequences alignments and the scores values, otherwise, only the latter are shown.</p>
        <p>the "Use pruning" uses the the pruning technique discussed in the article (<a href="http://bioinformatics.oxfordjournals.org/content/29/8/996" target="_blank">link</a>). The value corresponds to the lower bound set size. This option may reduce considerably the computation time if adjusted properly. Experimental results suggest that 10 is an appropriate value.</p>

        <br>
        <h3 id="step4">Step 4 - Submit</h3> <br>
        <p>This is the final step. The score values are displayed in a plot with indels/gaps as the y-axis and matches/score as the x-axis. It is possible to zoom in the plot by clicking it and drag. Double click to reset. If the option <b>traceback</b> was selected it is also possible to click on a data point in the plot and see the respective alignment below the plot. One can also compare all the alignments (button below plot). Identical (or similar) residues are given a color as it is described in <a href="http://www.bioinformatics.org/sms2/color_align_prop.html" target="_blank">http://www.bioinformatics.org/sms2/color_align_prop.html</a>. The coloring used is the default option of the website above. </p>
        <p>Links to download files with the score values and corresponding alignments are presented, so that the user can visualize the results again by uploading them in <a href="http://mosal.dei.uc.pt/align/" target="_blank">mosal.dei.uc.pt/align/</a> after selecting the option <b>Visualize results</b> in the top of page.</p>
        <p>For large sequences, the computation time can be very large. The user has the option of receiving the resulting files by e-mail, so that it can later visualize them.</p>

        <br>
        <div class="navbar">
            <?php include '../includes/footer.php'; ?>
        </div>
</body>
</html>

