<script>
if(window['cytoscape'] === undefined){
    var paths = {
	cytoscape: 'http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cytoscape.min.js',
    };

    console.log('Begin loading all JavaScript libs...');
    require.config({paths: paths});

    require(['cytoscape'], function(cytoscape){
	window['cytoscape'] = cytoscape;
	console.log('Finished loading jQuery and Cytoscape.js.');

	var event = document.createEvent("HTMLEvents");
	event.initEvent("load_cytoscape", true, false);
	window.dispatchEvent(event);
    });
}
</script>
