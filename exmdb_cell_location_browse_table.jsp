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

}
</style>



</head>

<% 
String searchname = "cytoplasm";  
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
	searchname = request.getParameter("searchname");
}

String table_sql = "select * from data_location_all where location = \""+ searchname + "\"";
// System.out.println(table_sql);

String search_table = dbhello.exmdb_location_browse_table(table_sql);

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
		         	 			<th rowspan="1" colspan="1" style="width:20px;">Gene name</th>	
		         	 			<th rowspan="1" colspan="1" style="width:40px;">Gene ID</th> 
		         	 			<th rowspan="1" colspan="1" style="width:40px;">Cell Location</th>
								<th rowspan="1" colspan="1" style="width:20px;">Tissue</th>
								<th rowspan="1" colspan="1" style="width:40px;">Method</th>
								<th rowspan="1" colspan="1" style="width:20px;">Pubmed id</th>
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
    	        '<option value="5">5</option>' +
    	        '<option value="10">10</option>' +
    	        '<option value="15">15</option>' +
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
    	     bAutoWidth: false,
//     	     "serverSide": true,         //开启服务端分页传输数据
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
<script src="assets/js/contact-form-script.js"></script>
<script src="assets/js/custom.js"></script>



</body>
</html>