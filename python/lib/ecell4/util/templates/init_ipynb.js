<script>
if(window['d3'] === undefined ||
   window['THREE'] === undefined ||
   window['Elegans'] === undefined){
    var paths = {
	d3: 'http://cdnjs.cloudflare.com/ajax/libs/d3/3.4.4/d3.min',
	THREE: 'http://cdnjs.cloudflare.com/ajax/libs/three.js/r66/three.min',
    };

    var shim = {
	'THREE': {exports: 'THREE'}
    };
    
    console.log('Begin loading all JavaScript libs...');
    require.config({paths: paths, shim: shim});

    require(['d3', 'THREE'], function(d3, THREE){
	window['d3'] = d3;
	window['THREE'] = THREE;
	console.log('Finished loading d3.js and Three.js.');

	// URL shown below (rawgit.com) should be replaced!
	var script = d3.select("head")
	    .append("script")
	    // .attr("src", 'https://rawgit.com/domitry/elegans/master/release/elegans.js')
	    .attr("src", 'https://raw.githubusercontent.com/domitry/elegans/microbial-sim/release/elegans.js')
	    .attr("async", true);

	script[0][0].onload = script[0][0].onreadystatechange = function(){
	    var event = document.createEvent("HTMLEvents");
	    event.initEvent("load_elegans", true, false);
	    window.dispatchEvent(event);
	    console.log('Finished loading Elegans.js ;)');
	};
    });
}
</script>
