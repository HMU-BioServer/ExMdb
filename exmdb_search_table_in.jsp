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

<style>
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



</style>

<!-- Datatable -->
<link rel="stylesheet" href="assets/css/style.min.css">
<link rel="stylesheet" type="text/css" href="Datatable/css/table.css">
<link rel="stylesheet" type="text/css" href="Datatable/css/buttons.dataTables.min2.css">
<link rel="stylesheet" href="Datatable/css/reset.css">
<link rel="stylesheet" href="Datatable/css/icon.css">

</head>

<% 
String searchname = "MIMAT0009451";  
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
	searchname = request.getParameter("searchname");
}

String organ = "lung";  
if(request.getParameter("organ")!=null&&!request.getParameter("organ").equals("")){
	organ = request.getParameter("organ");
}

String amount = "1";  
if(request.getParameter("amount")!=null&&!request.getParameter("amount").equals("")){
	amount = request.getParameter("amount");
	
}

String test_str = "";
test_str=dbhello.exmdb_search_trans(searchname);
test_str = test_str.substring(0,test_str.length() - 1);

String sql = "select * from organ_view where genename in ("+test_str+") and organ=\""+organ+"\" and fdr<="+amount;

// System.out.println("sql: "+sql);

String search_table = dbhello.search_table(sql);
String barset = dbhello.More_data(sql);
%>


<body>

<div id="Bar" style="height: 600px;width:50%;float:left"></div>
<div id="count_bar" style="height: 600px;width:50%;float:left""></div>


<div id="content-d" class="content light" style="width:100%;min-width:500px;float:left" >
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
		         	 			<th rowspan="1" colspan="1" style="width:60px;">Gene Name</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">P Value</th>
								<th rowspan="1" colspan="1" style="width:60px;">Data</th>
								<th rowspan="1" colspan="1" style="width:60px;">Organ</th>	
								<th rowspan="1" colspan="1" style="width:60px;">Detail</th>	
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



<script type="text/javascript" src="echart/echarts.min.js"></script>		
<script type="text/javascript">
var dom = document.getElementById("Bar");
var myChart = echarts.init(dom);
var app = {};

var option;



option = {
  legend: {textStyle:{fontFamily: 'Montserrat',fontSize:15}},
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15},},
  dataset: {
	  <%out.println(barset);%>
  },
  xAxis: {
    type: 'category',
      axisLabel: {
      fontFamily: 'Montserrat',
      rotate: 30
    }
  },
  yAxis: {},
  toolbox: {
	    feature: {
	      dataView: { show: true, readOnly: false },
	      restore: { show: true },
	      saveAsImage: { show: true }
	    }
	  },
  // Declare several bar series, each will be mapped
  // to a column of dataset.source by default.
  series: [
	       { type: 'bar',color: "#F5D043"}, 
	       { type: 'bar',color: "#114182"},
	       
	       ]
};




if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>

<script type="text/javascript">
var dom = document.getElementById("count_bar");
var myChart = echarts.init(dom);
var app = {};
var option;
option = {
  dataset: {
    source: [
      ['score', 'amount', 'product'],
      [1,1,'Bladder'],
      [1,1,'Oral '],
      [1,1,'Renal'],
      [1,1,'Tyriod'],
      [3,3,'Cholangiocarcinoma '],
      [3,3,'Liver'],
      [3,3,'Ovarian'],
      [4,4,'Bone'],
      [4,4,'Lung'],
      [4,4,'Prostate'],
      [5,5,'Glioma'],
      [5,5,'Pancreatic'],
      [6,6,'Breast'],
      [7,7,'Colorectal'],
      [7,7,'Gastric'],
      [10,10,'Esophageal'],
    ]
  },
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15},},
  grid: { 
      top: 50,
      bottom: 10,
      left: 0,
      right: 20,
      containLabel: true },
  xAxis: { name: ' ' },
  yAxis: { 
	  type: 'category',
      axisLabel: {fontFamily: 'Montserrat'}
  },
  visualMap: {
    orient: 'horizontal',
    left: 'center',
    top: 0,
    min: 0,
    max: 5,
    textStyle:{fontFamily: 'Montserrat',fontSize:15},
    text: ['','Data Counts'],
    // Map the score column to color
    dimension: 0,
    inRange: {
      color: ['#F5D043','#114182']
    }
  },
  series: [
    {
      type: 'bar',
          label: {
          show: true,
          precision: 2,
          position: 'right',
          valueAnimation: true,
          fontFamily: 'Montserrat'
        },
      encode: {
        // Map the "amount" column to X axis.
        x: 'amount',
        // Map the "product" column to Y axis
        y: 'product'
      }
    }
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>


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

</body>
</html>