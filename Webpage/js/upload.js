function executeUpload(){
	$(".alert").alert('close');
	var file = $("#data_file").get(0).files[0];  		// this field is required in html
	var reader = new FileReader();
	reader.readAsText(file, "UTF-8");
	reader.onload = function(evt){
		try{
			var D = parse_file_contents(evt.target.result);
			var header = evt.target.result.substring(0, evt.target.result.indexOf("\n"));
			var V = new Visualizer(D, 'results', header.indexOf("subscore")>=0, header.indexOf("indels")>=0);
			V.displayResults();
			$('#compare_all_aligns').click(function(){
					V.displayAllAligns();	
			});
		} catch(err){
			if(err === "parsing_file")
			handleError('Uploaded file <i>'+file.name+'</i> was not parsed correctly, please verify its contents (eg, extra spaces or newlines)!');
		}
	}
	reader.onerror = function (evt) {
		handleError('Uploaded file <i>'+file.name+'</i> could not be read, please try again!');
	}
}
function parse_file_contents(raw_data){
	// returns type of file (either values or alignments)
	var regexV = /^(matches|subscore)\t(indels|gaps)[\t ]*(\r?\n-?(0|[1-9][0-9]*)\t-?(0|[1-9][0-9]*)[\t ]*)+(\t| |\r?\n)*$/,
	regexA = /^(matches|subscore)\t(indels|gaps)[\t ]*(\r?\n-?(0|[1-9][0-9]*)\t-?(0|[1-9][0-9]*)[\t ]*\r?\n[A-Z\-\*]+[\t ]*\r?\n[A-Z\-\*]+[\t ]*)+(\t| |\r?\n)*$/;

	if( regexV.test(raw_data) )
		return get_values(raw_data);
	if( regexA.test(raw_data) )
		return get_aligns(raw_data);

	throw "parsing_file";
}
function get_values(raw_data){
return {values: raw_data.split("\n").slice(1,-1).map(parseInts)};
}
function get_aligns(raw_data){
var D = {values:[], alignments:[]}
	lines = remove_white_space(raw_data.split("\n"));

	for (var i = 1; i < lines.length-1; i+=3) {
		D.values.push(parseInts(lines[i]));			// push values
		D.alignments.push( [lines[i+1], lines[i+2]] );
	};
	return D;
}
function parseInts(str){
	arr = str.split("\t");
return arr.map(function(num){return parseInt(num);})
}
function handleError(error){
	var msg;
	switch(error){
		default:
		msg = error;
		break;
	}
	displayError('alerts_container', msg);
}
function remove_white_space(array){
	// removes white space from the end of the line
	new_array = [];
	for (var i = 0; i < array.length; i++) {
		new_array.push(array[i].replace(/[ \r\t]*$/, ''));
	};
	return new_array;
}