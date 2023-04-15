<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">


<link rel="stylesheet" href="tooltip_table/css/style.css" type="text/css">

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
.Detail_link{
	 transition: transform 0.25s ease;
	 background: #039a15;
	 padding: 12px;
	 border-radius: 5px;
	 border: solid 1px #039a15;
	 color:white !important;
	 cursor:pointer
 }
.Detail_link:hover {
	 background: white;
	 color:#039a15 !important;
	-webkit-transition: all 0.3s ease-in-out;
	-moz-transition: all 0.3s ease-in-out;
	-o-transition: all 0.3s ease-in-out;
	-ms-transition: all 0.3s ease-in-out;
	transition: all 0.3s ease-in-out;
}


div.material-table table tr td {
    text-align: center;
    line-height: 5px;
    height: 40px;
    font-size: 10px;
    color: rgba(0, 0, 0, 0.87);
    /* border-bottom: solid 1px #DDDDDD; */
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}
</style>



</head>

<%

String gseid="GSE142213";
if(request.getParameter("gseid")!=null&&!request.getParameter("gseid").equals("")){	
	gseid = request.getParameter("gseid");
}

String sql = "select * from tsne_cell_all where dataid =\"" +gseid + "\"";

String search_table = dbhello.exmdb_single_cell_table(sql);
// System.out.println(sql);

%>

<body>
<div id="content-d" class="content light" style="width:100%;min-width:500px;height:375px;float:left" >
		<div class="row scroll-pane">
		   <table id="datatable1" class="display nowrap" cellspacing="0" cellpadding="0" width="100%" >
		        <thead>
		         	 <tr>
		         	 			<th rowspan="1" colspan="1" style="width:60px;">Cell</th>	
		         	 			<th rowspan="1" colspan="1" style="width:60px;">angiogenesis</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">apoptosis</th>	
		         	 			<th rowspan="1" colspan="1" style="width:60px;">cellcycle</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">differentiation</th>	
		         	 			<th rowspan="1" colspan="1" style="width:60px;">dnadamage</th>	
		         	 			<th rowspan="1" colspan="1" style="width:60px;">dnarepair</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">emt</th>	
		         	 			<th rowspan="1" colspan="1" style="width:60px;">hypoxia</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">inflammation</th>
		         	 			<th rowspan="1" colspan="1" style="width:60px;">invasion</th>	
		         	 			<th rowspan="1" colspan="1" style="width:60px;">metastasis</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">proliferation</th>	
		         	 			<th rowspan="1" colspan="1" style="width:60px;">quiescence</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">stemness</th>
		            </tr> 
		        </thead>
		        <tbody>
		        <%out.println(search_table);%>
		        </tbody>
           	</table>
		</div>
</div>






       <script>
        var $tool = jQuery.noConflict();
		$tool(document).ready(function() {
			$tool('#datatable1').dataTable({
				"pagingType": "full_numbers",  //æ¾ç¤ºåé¡µææ
				"aoColumnDefs": [ { "bSortable": false, "aTargets": [ 0 ] }],
				"aaSorting": [],
	        "language": {
	            "lengthMenu": "Display _MENU_ records per page",
	            "info": "Showing page _PAGE_ of _PAGES_",
	            "infoEmpty": "No records available",
	            "infoFiltered": "(filtered from _MAX_ total records)"
	        },
			"lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],  //å®å¶é¿åº¦
			"scroller": false,
			dom: 'Bfrtip',
        buttons: [
            'copyHtml5',
            'excelHtml5',
            'csvHtml5',
            'pdfHtml5'
        ]
			});
		});
		
	  </script>

 <script src="http://www.jq22.com/jquery/1.11.1/jquery.min.js"></script>
 <script src="tooltip_table/js/tooltip_table.js"></script>

</body>
</html>