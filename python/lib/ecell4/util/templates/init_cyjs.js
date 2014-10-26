<script>
if(window['cytoscape'] === undefined){
    var paths = {
	cytoscape: 'http://cdnjs.cloudflare.com/ajax/libs/cytoscape/2.2.10/cytoscape'
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
