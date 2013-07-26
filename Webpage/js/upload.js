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
	var regexV = /^(matches|subscore)\t(indels|gaps)\n(-?(0|[1-9][0-9]*)\t-?(0|[1-9][0-9]*)\n)+$/,
		regexA = /^(matches|subscore)\t(indels|gaps)\n(-?(0|[1-9][0-9]*)\t-?(0|[1-9][0-9]*)\n[A-Z\-\*]+\n[A-Z\-\*]+\n)+$/;

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
	raw_data = raw_data.split("\n");

	for (var i = 1; i < raw_data.length-1; i+=3) {
		D.values.push(parseInts(raw_data[i]));			// push values
		D.alignments.push( [raw_data[i+1], raw_data[i+2]] );
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
		case "distinct_inputs":
			msg = "Please verify if the input is the same for both files (values and alignments).";
			break;
		default:
			msg = error;
			break;
	}
	displayError('alerts_container', msg);
}