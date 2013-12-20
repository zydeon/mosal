var ss_formdata = false;
var is_indels;
var is_traceback = false;
var is_subscore;
var wants_subscore;
var subscore_value;
var is_dp;			// pruning or not
var data = false;
var with_email = false;
var V;

$(document).ready(function() {
	$("#seq1").focus();
	$("#ibtn").button('toggle');
	$("#mbtn").button('toggle');
	$("#protein").button('toggle');

	$("#subscore_form").hide();
	$("#pruning_form").hide();
	id_dp = true;
	$('input[type=file]').bootstrapFileInput();
	$('#pruningValue').css('width',40);

	wants_subscore = false;
	// add events
	setEvents();
});

function setEvents(){
	$("#mbtn").click( function(){
		ss_formdata = false;
		wants_subscore = false;
		$("#subscore_form").hide();
		resetSubscoreSel($("#dna_select"));
		resetSubscoreSel($("#protein_select"));
		resetSubscoreUpl($("#subscore_upl"));
	} );
	$("#sbtn").click( function(){
		wants_subscore = true;
		$("#subscore_form").show();
	} );
	$("#pruning_chkbox").click( function(){
		if(this.checked){
			is_dp = false;
			$("#pruning_form").show();
			$("#pruningValue").val(10);
			$("#pruning").val(10);			
		}
		else{
			is_dp = true;
			$("#pruning_form").hide();
		}
	} );
	$("#pruning").change( function(){
		$("#pruningValue").val($(this).val());
	});
	$("#pruningValue").change( function(){
		if($(this).val() > 25)
			$(this).val(25);
		else if($(this).val() < 1)
			$(this).val(1);
		$("#pruning").val($(this).val());
	});
	$("#email_chkbox").click(function(){
		if(this.checked){
			$("#email_form").html( emailForm() );
			with_email = true;
			$("#email").focus();
		}
		else{
			$("#email_form").html("");
			with_email = false;
		}
	});
	$("#subscore_upl").click(function(){
		resetSubscoreSel($("#subscore_select"));
		if (window.FormData) ss_formdata = new FormData();
		else alert("Error: It's not possible to upload the substitution score file in this browser (does not support FormData)!");
	});
	$("#subscore_select").change(function(){
		ss_formdata = false;
		resetSubscoreUpl($("#subscore_upl"));
	});
	$("#dna").click(function(){
		$("#dna_select").show();
		$("#protein_select").hide();
	});
	$("#protein").click(function(){
		$("#dna_select").hide();
		$("#protein_select").show();
	});	
	$("#sample1").click(function(){
		if(isBtnActive("dna"))
			$("#seq1").val("> Sequence example 1\nATGAACAATCAAGCATACGGTGTTACACCTCCTATATCTGTTGCCAACTCTACTCCTAAG\nGAAAACGAACTTAATGATCTGTTAATAAAGGAATTGAAATCAAGAGGCTCGTTTGAGAGT\nGAAACTGCCACGAAGAAAAGAGTGGAGGTGTTAAATATTCTTCAAAGTATGACTGAAGAA\nTTTGTCTACAAAGTTTCTATAAAGAAAAACATATCGGAAGGAATGGCAAGGGATGTGGGT\nGGAAAAATATTTACATTTGGTTCCTATAGGTTAGGGGTATATGGACCAGGTTCAGATATT\nGACACTCTAGTTGTTGTTCC");
		else
			$("#seq1").val("> Sequence example 1\nMATGANATPLDFPSKKRKRSRWNQDTMEQKTVIPGMPTVIPPGLTREQERAYIVQLQIED\nLTRKLRTGDLGIPPNPEDRSPSPEPIYNSEGKRLNTREFRTRKKLEEERHNLITEMVALN\nPDFKPPADYKPPATRVSDKVMIPQDEYPEINFVGLLIGPRGNTLKNIEKECNAKIMIRGK\nGSVKEGKVGRKDGQMLPGEDEPLHALVTANTMENVKKAVEQIRNILKQGIETPEDQNDLR\nKMQLRELARLNGTLREDDNRILRPWQSSETRSI");
	});
	$("#sample2").click(function(){
		if(isBtnActive("dna"))
			$("#seq2").val("> Sequence example 2\nATGAACACGAAGACATACGGAGTAACTGAGCCTATATCAACAAATGGTCCAACACCGAAG\nGAGAATATTCTCAATGATGCACTCATACAAGAATTGAAAAATAGAGGTTCGTTTGAAAGC\nGAACAAGCAACTAAAAAAAGAGTTGAAGTGTTGACATTATTTCAACGATTAGTTCAAGAA\nTTTGTGTACACAGTATCAAAGAGCAAAAACATGTCCGATTCAATGGCCCAAGATGCCGGA\nGGAAAGGTATTTACTTTTGGATCCTATAGGTTGGGAGTCTATGGTCCAGGATCTGATATT\nGACACATTGGTTGTTGTGCC");
		else
			$("#seq2").val("> Sequence example 2\nMATKTLARPDEPGVHARDPRSLHPQDETSKIVGNSASIPSDNNEPPTNSNEFESAEHSSK\nKSLVQIATVMASLCACVFLAALEVTIVSTALPTIAAHFTSDSGYTWIGTSFVLAHTASTP\nSWGKISDIWGRKPILLIANVIFFAGSLLCALVDDLAIFIAGRAIQGLGAAVWAVASGVGP\nILGGAFTVRLSWRWCFWINLPITVAVFFLLVLTLRLPSPNTPVWAGLKAIDWPGSFLIVG\nGTLMLLLGLYLGGVYEPWNSATVVCLIVFGIITALLFVWNEWKLAEYPVIPVHLFKTWSS");
	});	
}
function resetSubscoreUpl(element){
	element.replaceWith( element = element.clone( true ) );
	$('#subscore_form > span.file-input-name').remove();
}
function resetSubscoreSel(element){
	element.val('-');
}
function startAlignment(){
	var timestamp = getTimestamp();

	// reset alerts
	$(".alert").alert('close');

	// validate sequences
	try{
		validateSeq($("#seq1").val())
	} catch(err){
		displayError('alerts_seq1', err);
		return;
	}
	try{
		validateSeq($("#seq2").val())
	} catch(err){
		displayError('alerts_seq2', err);
		return;
	}

	// upload substitution score file
	if(wants_subscore){
		var files = $("#subscore_upl").get(0).files;			// field is required (length is always > 0)
		if(ss_formdata && files.length > 0 ){
			var f = files[0];
			var reader = new FileReader();
			reader.readAsText(f, "UTF-8");
			reader.onload = function(evt){
				// validate file first
				if(validSubscoreFile(evt.target.result)){
					ss_formdata.append("subscore_file", f);
					ss_formdata.append("timestamp", timestamp);

					$.ajax({
						url: "upload_subscore.php",
						type: "POST",
						data: ss_formdata,
						processData: false,
						contentType: false,
						success: function (response) {
							var response = JSON.parse(response);
							if(response.status == 'error'){
								displayError('alerts_container','Substitution score file <i>'+f.name+'</i> upload failed.');
							}
							else if(response.status == 'success') {
								subscore_value = "upload";
								executeAlign(timestamp);
							}
						}
					});
				}
				else{
					displayError('alerts_subscore','Please verify syntax of the substitution score file <i>'+f.name+'</i>.');
				}
			}
			reader.onerror = function (evt) {
				displayError('alerts_subscore','Uploaded substitution score file <i>'+f.name+'</i> could not be read, please try again!');
			}
		}
		else{
			var is_dna = isBtnActive("dna");
			var sel_value = $(is_dna ? "#dna_select" : "#protein_select").val();
			if(sel_value != "-"){
				subscore_value = sel_value;
				executeAlign(timestamp);
			}
			else
				displayError('alerts_subscore','Please make sure the substitution matrix chosen is correct!');
		}
	}
	else{
		subscore_value = false;
		executeAlign(timestamp);
	}
}
function executeAlign(timestamp){
	// reset graphs and set loading
	$("#load").html("<img src='../img/loading.gif' alt='loading...' />");
	$("#results").html("");

	// compute alignment
	is_indels    = isBtnActive("ibtn");
	is_subscore  = isBtnActive("sbtn");
	is_traceback = $("#tb_chkbox").is(':checked');
	var bound_number = 0;
	if(is_dp) bound_number = $("#pruningValue").val();

	var email = with_email ? $("#email").val() : "",
		email_subject = with_email ? ($("#subject").val()!="" ? $("#subject").val() : "MOSAL sequence alignment solutions") : "";

	$.post(
		"align.php",
		{	
			seq1 		: clean_sequence($("#seq1").val()),
			seq2 		: clean_sequence($("#seq2").val()),
			timestamp 	: timestamp,
			is_indels 	: is_indels,
			bound_number: bound_number,
			is_traceback: is_traceback,
			subscore    : subscore_value,
			email 		: email,
			email_subject 		: email_subject,

		},
		function(response){
			var response = JSON.parse(response);
			$("#load").html("");

			if(response.status=='success'){
				V = new Visualizer(response.data, 'results', is_subscore, is_indels)
				V.displayResults();
				$('#compare_all_aligns').click(function(){
						V.displayAllAligns();	
				});
				displayDLButtons();
				if(response.debug)
					console.log("Debug info:\n> "+response.debug);
			}
			else if(response.status=='error'){
				console.log(response.info);
				console.log("Debug info:\n> "+response.debug);
			}
		}
	);	
}
function validateSeq(value){
	var regex;
	if(value.indexOf(">") != -1){			// fasta format
		regex = /^(\t| |\r?\n)*>.*(\r?\n[A-Z\*]+[\t ]*)+(\t| |\r?\n)*$/;
		if(!regex.test(value)) throw "Please make sure the sequence below is either accordingly FASTA format or just a string of characters (only values between A-Z and the * character are allowed)"
	}
	else{
		regex = /^(\t| |\r?\n)*[A-Z\*]+(\t| |\r?\n)*$/;
		if(!regex.test(value)) throw 'Please check the sequence below for illegal characters: only values between A-Z and the * character are allowed'
	}
}
function displayDLButtons(){
	var html = 	'<a id="downloadValues" download="values.tsv" class="btn btn-success"><i class="icon-white icon-circle-arrow-down"></i> Download all solutions</a>' +
				'<br>';
	if(is_traceback)
		html += '<br>' +
				'<a id="downloadAligns" download="alignments.txt" class="btn btn-success"><i class="icon-white icon-circle-arrow-down"></i> Download all alignments</a>';

	$("#results").append(html);
	$("#downloadValues").click( function(){
		var s1 = is_subscore ? "subscore" : "matches",
			s2 = is_indels ? "indels" : "gaps";
		$(this).attr('href','data:text/tab-separated-values;base64,'+btoa(V.dataValuesToCSV(s1+"\t"+s2)));
	});
	$("#downloadAligns").click( function(){
		var s1 = is_subscore ? "subscore" : "matches",
			s2 = is_indels ? "indels" : "gaps";
		$(this).attr('href','data:text/plain;base64,'+btoa(V.dataAlignsToTxt(s1+"\t"+s2)));
	});		
}
function isBtnActive(id){
	return $("#"+id).hasClass('active');
}
function getTimestamp(){
	return new Date().getTime();
}
function emailForm(){
	return 	'<table>' +
				'<tr>' +
					'<td class="right_align">E-mail address</td>' +
					'<td><input type="email" id="email" required></td>' +
				'</tr>' +
				'<tr>' +
					'<td class="right_align">Subject</td>' +
					'<td><input type="text" id="subject"></td>' +
				'</tr>' +
			'</table>';
}
function validSubscoreFile(raw_data){
	var regex = /^(#.*\r?\n)*([\t ]*[A-Z\*])+[\t ]*(\r?\n[A-Z\*]([\t ]*-?(0|[1-9][0-9]*))+[\t ]*)+(\t| |\r?\n)*$/;
	return regex.test(raw_data);
}
function clean_sequence(seq){	
	seq = seq.replace(/^(\t| | \r?\n)*/,'');	// beginning of input
	seq = seq.replace(/[\t ]*\r?\n/g,'\r\n');	// end of lines
	return seq.replace(/(\t| | \r?\n)*$/,'');	// end of input
}
