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
            top:100px;
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
        <p>This webpage describes the process of using <a href="http://mosal.dei.uc.pt/align">MOSAL</a> to get multiobjective sequence alignments. In order to use the tool, one must follow a multiple step process that is explained below.</p>
        <br>
        <h3 id="step1">Step 1 - Sequences input</h3> <br>
        <p>This is the step where the user provide the target sequences. The format of the sequences can either be the <a href="http://en.wikipedia.org/wiki/FASTA_format" target="_blank">FASTA</a> format or just a string of characters. In both formats only uppercase letters are allowed together with the * symbol. </p>
        <p>The choice between <b>Protein</b> and <b>DNA</b> influences the choice of the substitution score matrix (step 2). There are also buttons below each sequence input field that introduce an example of a FASTA formatted sequence.</p>

        <br>
        <h3 id="step2">Step 2 - Score function</h3> <br>
        <p>In this step the user must choose 2 objectives:</p>
        <ul>
            <li>Minimize either <b>Indels</b> or <b>Gaps</b></li>
            <li>Maximize either the number of <b>Matches</b> or the <b>Substitution score</b> </li>
        </ul>
        <p>The user can choose one of the default substitution score matrixes (dependent on choice between Protein or DNA) or upload one of its own as long as it follows the standards (see <a href="http://en.wikipedia.org/wiki/Substitution_matrix" target="_blank">link</a>).</p>

        <br>
        <h3>Step 3 - Alignment options</h3> <br>
        <p>The first option in this section, when selected, tells the tool to also compute the sequences alignments besides computing the values of the objectives (default behaviour). Naturally, the time it will take to finish the execution is larger, but one can compare the alignments by clicking on the solutions in the chart that is displayed as output.</p>
        <p>The second option may reduce the computation time if adjusted properly: too few bounds may not increase the speed, but too many bounds may introduce more overhead during pruning.</p>

        <br>
        <h3>Step 4 - Submit</h3> <br>
        <p>This is the final step. The computation values of the objectives are displayed in a plot with indels/gaps as the y-axis and matches/score as the x-axis. It is possible to zoom in the plot by clicking it and drag. Double click to reset. If the option <b>traceback</b> was selected it is also possible to click on a data point in the plot and see the respective alignment below the plot. </p>
        <p>Links to download files with the computed values and alignments are presented, so that the user can visualize the results again by uploading them in <a href="http://mosal.dei.uc.pt/align/">mosal.dei.uc.pt/align/</a> after selecting the option <b>Visualize results</b> in the top of page.</p>
        <p>For large sequences, the time of computation can be very large, so the the user has the option of receiving the result files mentioned by e-mail, so that it can later visuzalize them.</p>

        <br>
        <div class="navbar">
            <?php include '../includes/footer.php'; ?>
        </div>
</body>
</html>

