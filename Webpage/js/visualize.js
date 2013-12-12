function Visualizer(data, targetDivID, is_subscore, is_indels){
    this.data = data;

    this.is_traceback = !!data.alignments;
    this.is_subscore  = is_subscore;
    this.is_indels    = is_indels;
    this.targetDivID  = targetDivID;

    this.displayResults = function(){
        $("#"+targetDivID).html( this.resultsHTML() );

        this.setEvents();

        scrollTo(this.targetDivID);
        this.displayChart(this.data.values);
    }
    this.setEvents = function(){
        if(this.is_traceback){
            $('#chartdiv').bind('jqplotDataClick',
                    function (ev, seriesIndex, pointIndex, data) {                
                        $("#aligndiv").html( createAlignDiv(pointIndex) );
                        scrollTo(targetDivID);
                    }
                )
            .bind('jqplotHighlighterHighlight', function () {
                $('.jqplot-event-canvas').css( 'cursor', 'pointer' );
            }).bind('jqplotHighlighterUnhighlight', function () {
                $('.jqplot-event-canvas').css( 'cursor', 'crosshair' );
            });     
        }
    }
    function createAlignDiv(index){
        var x  = data.values[index][0], 
            y  = data.values[index][1],
            a1 = data.alignments[index][0],
            a2 = data.alignments[index][1];

        var b = "";
        for (var i = 0; i < a1.length; i++)
            b += a1[i]==a2[i] ? "|" : "&nbsp;";

        return $("<div class='align'></div>")
                .html(  "<br><pre class='align'>" +
                            color (a1) + "<br>" + 
                            b + "<br>" +
                            color (a2) + "<br>" +
                        "</pre>" +
                        "("+x + ", " + y + ")<br>"
                    );
    }
    function color(sequence){
        color_ = getColorsMap();
        return sequence.split("").map(function(c){
            return "<span style='color:"+color_[c]+"'>" + c + "</span>";
        }).join("");
    }
    function getColorsMap(){
        M   = {} // maps letters to colors
        AB  = "-ABCDEFGHIJKLMNOPQRSTUVWXYZ*".split("");
        inc = Math.floor(256/(Math.pow(28, 1/3)));   // get equally (maximum) spread colors
        r = 0;
        g = 0;
        b = 0;
        i = 0;
        for (var r = 0; r < 256 && i < AB.length; r+=inc) {
            for (var g = 0; g < 256 && i < AB.length; g+=inc) {
                for (var b = 0; b < 256 && i < AB.length; b+=inc) {
                    console.log(AB[i]+": "+toHexRepr(r,g,b))
                    M[ AB[i] ] = toHexRepr(r,g,b);
                    i++;   // guarantee that does not exceed number of colors
                }
            }
        }
        return M;
    }
    function toHexRepr(r,g,b){
        r_ = (r > 15 ? "" : "0") + r.toString(16);
        g_ = (g > 15 ? "" : "0") + g.toString(16);
        b_ = (b > 15 ? "" : "0") + b.toString(16);
        return "#" + r_ + g_ + b_;
    }
    this.resultsHTML = function(){
        var contents =  '<h2>Results</h2>' + 
                        '<br>' +
                        '<div id="chartdiv"style="height:400px;width:600px;"></div>' +
                        '<br>' +
                        '<div id="aligndiv"></div>' +
                        '<br>';

        return contents;
    }
    this.displayChart = function(data){
        var xlabel = this.is_subscore ? 'Substitution score' : '# Matches',
            ylabel = this.is_indels   ? '# Indels'           : '# Gaps';

        // ---------- just in case data has only one point renderer does not act well
        var xmin = xmax = ymin = ymax = undefined;
        if(data.length===1){
            xmin = data[0][0]-1;
            xmax = data[0][0]+1;
            ymin = data[0][1]-1;
            ymax = data[0][1]+1;
        }
        // ---------------

        $.jqplot('chartdiv',  [data], {
            title:data.length+' non-dominated alignments found',
            seriesColors:['#006dcc'],
            axes:{
                xaxis:{
                    min: xmin,
                    max: xmax,
                    label: xlabel,
                    labelRenderer: $.jqplot.CanvasAxisLabelRenderer,
                    labelOptions: {
                        fontFamily: 'Helvetica Neue',
                        fontSize: '12pt'
                    }
                },
                yaxis:{
                    min: ymin,
                    max: ymax,
                    label: ylabel,
                    labelRenderer: $.jqplot.CanvasAxisLabelRenderer,
                    labelOptions: {
                        fontFamily: 'Helvetica Neue',
                        fontSize: '12pt'
                    }
                }
            },
            cursor:{
                show: true,
                zoom: true,
                showTooltip: false,
            },
            highlighter:{
                show: true,
                sizeAdjust: 15,
                formatString: '(%s , %s)',
            },
        });
    }
    this.downloadValues = function(){
        var s1 = this.is_subscore ? "subscore" : "matches",
            s2 = this.is_indels ? "indels" : "gaps";
        $(this).attr('href','data:text/plain;base64,'+btoa(dataAlignsToTxt(s1+"\t"+s2, this.data)));        
    }
    this.dataValuesToCSV = function(header){
        var text = header+"\n";
        for (var i = 0; i < this.data.values.length; i++) {
            text += this.data.values[i].join("\t") + "\n";
        };
        return text;
    }
    this.dataAlignsToTxt = function(header){
        var text = header+"\n";
        for (var i = 0; i < data.alignments.length; i++) {
            text += this.data.values[i].join("\t") + "\n" +
                    this.data.alignments[i].join("\n") + "\n";
        };
        return text;
    }

}
function displayError(divid, msg){
    $("#"+divid).append(
            '<div class="alert alert-error fade in">' +
                '<button type="button" class="close" data-dismiss="alert">&times;</button>' +
                '<strong>Error!</strong> ' + msg +
            '</div>'
        );
    scrollTo(divid);
}
function scrollTo(id){
    $('html,body').animate( {scrollTop: $("#"+id).offset().top },
                            'slow' );
}