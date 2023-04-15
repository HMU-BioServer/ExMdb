<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">



<link rel="stylesheet" href="assets/css/bootstrap.min.css">
<link rel="stylesheet" href="assets/css/animate.min.css">
<link rel="stylesheet" href="assets/fonts/flaticon.css">
<link rel="stylesheet" href="assets/css/boxicons.min.css">
<link rel="stylesheet" href="assets/css/owl.carousel.min.css">
<link rel="stylesheet" href="assets/css/owl.theme.default.min.css">
<link rel="stylesheet" href="assets/css/magnific-popup.css">
<link rel="stylesheet" href="assets/css/nice-select.min.css">
<link rel="stylesheet" href="assets/css/meanmenu.css">
<link rel="stylesheet" href="assets/css/style.css">
<link rel="stylesheet" href="assets/css/responsive.css">
<link rel="icon" type="image/png" href="assets/images/favicon.png">


<!-- Datatable -->
<link rel="stylesheet" href="assets/css/style.min.css">
<link rel="stylesheet" type="text/css" href="Datatable/css/table.css">
<link rel="stylesheet" type="text/css" href="Datatable/css/buttons.dataTables.min2.css">
<link rel="stylesheet" href="Datatable/css/reset.css">
<link rel="stylesheet" href="Datatable/css/icon.css">

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

 .evenflow_scale{
	-webkit-transition: all 0.3s ease-in-out;
	-moz-transition: all 0.3s ease-in-out;
	-o-transition: all 0.3s ease-in-out;
	-ms-transition: all 0.3s ease-in-out;
	transition: all 0.3s ease-in-out;
 }
.evenflow_scale:hover {
    transform: scale(1.3,1.3);
}

body::-webkit-scrollbar-track
{
/* 	-webkit-box-shadow: inset 0 0 6px rgba(0,0,0,0.1); */
	background-color: white;
	border-radius: 10px;
}

body::-webkit-scrollbar
{
	width: 5px;
	background-color: white;
}

body::-webkit-scrollbar-thumb
{
	border-radius: 10px;
	width: 5px;
	height: 5px;
	background-color: white;
}


.Detail_link {
    transition: transform 0.25s ease;
    background: #f5bf29;
    padding: 12px;
    border-radius: 5px;
    border: solid 1px #f5bf29;
    color: white !important;
    cursor: pointer;
}

.Detail_link:hover {
	 background: white;
	 color:#000000 !important;
	-webkit-transition: all 0.3s ease-in-out;
	-moz-transition: all 0.3s ease-in-out;
	-o-transition: all 0.3s ease-in-out;
	-ms-transition: all 0.3s ease-in-out;
	transition: all 0.3s ease-in-out;
}


.Detail_link_2 {
    transition: transform 0.25s ease;
    background: #2984f5;
    padding: 12px;
    border-radius: 5px;
    border: solid 1px #2984f5;
    color: white !important;
    cursor: pointer;
}

.Detail_link_2:hover {
	 background: white;
	 color:#000000 !important;
	-webkit-transition: all 0.3s ease-in-out;
	-moz-transition: all 0.3s ease-in-out;
	-o-transition: all 0.3s ease-in-out;
	-ms-transition: all 0.3s ease-in-out;
	transition: all 0.3s ease-in-out;
}
</style>



</head>

<% 
String searchname = "MALAT1";  


if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
	searchname = request.getParameter("searchname");
}

String search_table = dbhello.getSC_datainfo(searchname);

%>


<body>

<div class="preloader">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>

<div id="content-d" class="content light" style="width:100%;min-width:500px;height:375px;float:left" >
		<div class="row scroll-pane">
		  <div id="admin" class="col s12">
		    <div class="card material-table">
		      <div class="table-header" id="buttons">
		        <div class="actions">
		          <a class="buttons-toggle waves-effect btn-flat nopadding"   href="javascript:;" οnclick="redirect($(this))" val="#"><i class="material-icons">folder_shared</i></a>
		          <a class="search-toggle waves-effect btn-flat nopadding" href="javascript:;" οnclick="redirect($(this))" val="#"><i class="material-icons">search</i></a>
		        </div>
		      </div>
		      <table id="datatable1">
		        <thead>
		         	 <tr>
		         	 			
		         	 			<th rowspan="1" colspan="1" style="width:120px;">Data ID</th>
		         	 			<th rowspan="1" colspan="1" style="width:40px;">Disease</th> 
		         	 			<th rowspan="1" colspan="1" style="width:40px;">Sample</th> 
		         	 			<th rowspan="1" colspan="1" style="width:40px;">organ</th> 	
		         	 			<th rowspan="1" colspan="1" style="width:40px;">Species</th> 
		         	 			<th rowspan="1" colspan="1" style="width:80px;">Cite</th>
								<th rowspan="1" colspan="1" style="width:40px;">Pubmed id</th>
								<th rowspan="1" colspan="1" style="width:40px;">Detail</th>
		            </tr> 
		        </thead>
		        <tbody>
		        <%out.println(search_table);%>
		        </tbody>
		      </table>
		    </div>
		  </div>
		</div>
</div>




<script src="Datatable/js/jquery.min.js"></script>
<script type="text/javascript" language="javascript" src="Datatable/js/jquery.dataTables.js"></script>
<script src="Datatable/js/materialize.min.js"></script>
<script src="Datatable/js/table.js"></script>
<script type="text/javascript" language="javascript" src="Datatable/js/dataTables.buttons.min.js"></script>
<script type="text/javascript" language="javascript" src="Datatable/js/jszip.min.js"></script>
<script type="text/javascript" language="javascript" src="Datatable/js/pdfmake.min.js"></script>
<script type="text/javascript" language="javascript" src="Datatable/js/vfs_fonts.js"></script>
<script type="text/javascript" language="javascript" src="Datatable/js/buttons.html5.min.js"></script>

   <script>
      $(document).ready(function() {
    	  var table=$('#datatable1').dataTable({
    	    "oLanguage": {
    	      "sStripClasses": "",
    	      "sSearch": "",
    	      "sSearchPlaceholder": "Enter Keywords Here",
    	      "sInfo": "_START_ -_END_ of _TOTAL_",
    	      "sLengthMenu": '<span>Rows per page:</span><select class="browser-default">' +
    	        '<option value="10">10</option>' +
    	        '</select></div>'
    	    },
	        buttons: [
	            'copyHtml5',
	            'excelHtml5',
	            'csvHtml5',
	            'pdfHtml5'
	        ],
	        "aaSorting": [],
    	    "Processing":true,
    	     bAutoWidth: false
    	  }); 
    	  //高亮鼠标所在行和列
      });
  </script>


<script src="assets/js/owl.carousel.min.js"></script>
<script src="assets/js/jquery.magnific-popup.min.js"></script>
<script src="assets/js/jquery.nice-select.min.js"></script>
<script src="assets/js/wow.min.js"></script>
<script src="assets/js/meanmenu.js"></script>
<script src="assets/js/jquery.ajaxchimp.min.js"></script>
<script src="assets/js/form-validator.min.js"></script>
<script src="assets/js/custom.js"></script>


</body>
</html>