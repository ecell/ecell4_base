<script>
if(window['d3'] === undefined ||
   window['THREE'] === undefined ||
   window['Elegans'] === undefined){
    var paths = {
	d3: 'http://cdnjs.cloudflare.com/ajax/libs/d3/3.4.4/d3.min',
	THREE: 'http://cdnjs.cloudflare.com/ajax/libs/three.js/r66/three.min',
	Elegans: 'https://rawgit.com/domitry/elegans/master/release/elegans'
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
	require(['Elegans'], function(Elegans){
	    window['Elegans'] = Elegans;
	    console.log('Finished loading Elegans.js');
	    for(var key in paths){
		d3.select('head')
		    .append('script')
		    .attr('type', 'text/javascript')
		    .attr('src', paths[key] + '.js');
	    }
	    console.log('Finished loading ;)');
	});
    });
}
</script>
