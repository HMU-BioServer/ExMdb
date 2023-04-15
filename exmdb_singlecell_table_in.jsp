<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">


<link rel="stylesheet" href="tooltip_table/css/style.css" type="text/css">
<link rel="stylesheet" href="assets/css/style.css">
<!-- Datatables -->
  <script src="js_search_quick/hm.js"></script>
  <script type="text/javascript" src="js_search_quick/modernizr.custom.79639.js"></script>
  <link rel="stylesheet" type="text/css" href="js_search_quick/jquery.dataTables.min.css">

  <link rel="stylesheet" type="text/css" href="js_search_quick/buttons.dataTables.min2.css">
  <script type="text/javascript" language="javascript" src="js_search_quick/jquery.js"></script>
  <script type="text/javascript" language="javascript" src="js_search_quick/jquery.dataTables.js"></script>
  <script type="text/javascript" language="javascript" src="js_search_quick/dataTables.buttons.min.js"></script>
  <script type="text/javascript" language="javascript" src="js_search_quick/jszip.min.js"></script>
  <script type="text/javascript" language="javascript" src="js_search_quick/pdfmake.min.js"></script>
  <script type="text/javascript" language="javascript" src="js_search_quick/vfs_fonts.js"></script>
  <script type="text/javascript" language="javascript" src="js_search_quick/buttons.html5.min.js"></script>




<style>
.top_head {
    position: relative;
}

.top_head .tooltiptext {
    visibility: hidden;
    transition: transform 0.25s ease;
    margin: -5px 0px 0 5px;
    background-color: black;
    color: #fff;
    text-align: center;
    border-radius: 6px;
    padding: 5px 10px;

    /* 定位 */
    position: absolute;
    z-index: 1;
}

.top_head:hover .tooltiptext {
    visibility: visible;
    -webkit-transition: all 0.3s ease-in-out;
	-moz-transition: all 0.3s ease-in-out;
	-o-transition: all 0.3s ease-in-out;
	-ms-transition: all 0.3s ease-in-out;
	transition: all 0.3s ease-in-out;
}
</style>



</head>

<%

String gseid="GSE83142";
if(request.getParameter("gseid")!=null&&!request.getParameter("gseid").equals("")){	
	gseid = request.getParameter("gseid");
}

String sql = "select * from tsne_cell_all where dataid =\"" +gseid + "\"";

String search_table = dbhello.exmdb_single_cell_table(sql);
// System.out.println(sql);

%>

<body>

		   <table id="datatable1" class="display nowrap" cellspacing="0" cellpadding="0" width="100%" style="text-align: center;border-top: 1px solid #d2d2d2;border-bottom: 1px solid rgb(210, 210, 210);">
		        <thead>
		         	 <tr>
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:100px;"><img src="assets/images/icon/cell.png"> <span class="tooltiptext">Cell [Cell name]</span></th>	
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:80px;"><img src="assets/images/icon/celltype.png"> <span class="tooltiptext">Cell Type</span></th>	
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:80px;"><img src="assets/images/icon/sample.png"> <span class="tooltiptext">Sample</span></th>	
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:80px;"><img src="assets/images/icon/source.png"> <span class="tooltiptext">Source</span></th>	
		         	 			
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/angiogenesis.png"><span class="tooltiptext">Angiogenesis</span> </th>	
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/apoptosis.png"><span class="tooltiptext">Apoptosis</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/cellcycle.png"><span class="tooltiptext">Cell Cycle</span> </th>	
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/differentiation.png"><span class="tooltiptext">Differentiation</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/dnadamage.png"><span class="tooltiptext">DNA Damage</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/dnarepair.png"><span class="tooltiptext">DNA Repair</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/EMT.png"><span class="tooltiptext">EMT</span> </th>	
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/hypoxia.png"><span class="tooltiptext">Hypoxia</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/inflammation.png"><span class="tooltiptext">Inflammation</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/invasion.png"><span class="tooltiptext">Invasion</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/meta.png"><span class="tooltiptext">Metastasis</span> </th>	
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/proliferation.png"><span class="tooltiptext" style="margin-left: -195%;">Proliferation</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/quiescence.png"><span class="tooltiptext" style="margin-left: -195%;">Quiescence</span> </th> 
		         	 			<th class="top_head" rowspan="1" colspan="1" style="width:40px;padding:0"><img src="assets/images/icon/stemness_title.png"><span class="tooltiptext" style="margin-left: -195%;">Stemness</span> </th> 
		            </tr> 
		        </thead>
		        <tbody>
		        <%out.println(search_table);%>
		        </tbody>
           	</table>






       <script>
        var $tool = jQuery.noConflict();
		$tool(document).ready(function() {
			$tool('#datatable1').dataTable({
				"pagingType": "full_numbers",  //æ¾ç¤ºåé¡µææ
				"aoColumnDefs": [ { "bSortable": false, "aTargets": [ 0,1,2,3 ] }],
				"aaSorting": [],
	        "language": {
	            "lengthMenu": "Display _MENU_ records per page",
	            "info": "Showing page _PAGE_ of _PAGES_",
	            "infoEmpty": "No records available",
	            "infoFiltered": "(filtered from _MAX_ total records)"
	        },
			"lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],  //å®å¶é¿åº¦
			"scroller": false,
			dom: 'Brtip',
        buttons: [
            'copyHtml5',
            'excelHtml5',
            'csvHtml5',
            'pdfHtml5'
        ]
			});
		});
		
	  </script>

 <script src="js/1.11.1/jquery.min.js"></script>
 <script src="tooltip_table/js/tooltip_table.js"></script>

</body>
</html>