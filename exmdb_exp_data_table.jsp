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

</style>



</head>

<% 
String searchname = "RP11-838N2.4";  
String expvalStr = "";
String expvalStr_sub = "";

if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
	searchname = request.getParameter("searchname");
}

if(request.getParameter("expval")!=null&&!request.getParameter("expval").equals("")){
	expvalStr = request.getParameter("expval");
}

for (String retval: expvalStr.split(",")){
	String tempStr =  "method like '%"+retval+"%'";
	tempStr = tempStr + " or ";
	expvalStr_sub = expvalStr_sub + tempStr;
}

expvalStr_sub = expvalStr_sub.substring(0, expvalStr_sub.length() -3);	
expvalStr_sub="("+expvalStr_sub+")";



String table_sql = "select * from experimental_data where (genename = \""+searchname+"\" or cancername =\"" + searchname + "\") and "+ expvalStr_sub;
// System.out.println(table_sql);

String search_table = dbhello.exmdb_exp_table(table_sql);

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
		         	 			<th rowspan="1" colspan="1" style="width:40px;">Gene name</th>	
		         	 			<th rowspan="1" colspan="1" style="width:40px;">Cancer name</th> 
		         	 			<th rowspan="1" colspan="1" style="width:120px;">Method</th>
								<th rowspan="1" colspan="1" style="width:60px;">Pubmed id</th>
								<th rowspan="1" colspan="1" style="width:20px;">Year</th>
								<th rowspan="1" colspan="1" style="width:20px;">Detail</th>
		            </tr> 
		        </thead>
		        <tbody>
		        <%out.println(search_table);%>
		        </tbody>
		      </table>
		    </div>
		  </div>
		</div>
		
		
		<div id="container_bar" style="height: 400px"></div>
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

<script type="text/javascript" src="echart/echarts.min.js"></script>
<script type="text/javascript">
var dom = document.getElementById("container_bar");
var myChart = echarts.init(dom);
var app = {};

var option;

option = {
		  title: {text: 'Data Statistics',
			   textStyle: {
				   fontSize: 15,
				   fontFamily: 'Montserrat'
			   }},
  dataset: {
    source: [
      ['score', 'amount', 'product'],
      [527,527,'Others .. [Count < 20]'],
      [20,20,'Osteosarcoma'],
      [24,24,'Cervical Cancer'],
      [26,26,'papillary thyroid cancer'],
      [28,28,'Liver Cancer'],
      [28,28,'Osteoarthritis'],
      [28,28,'Liver Fibrosis'],
      [30,30,'Esophageal Cancer'],
      [36,36,'Alzheimer Disease'],
      [45,45,'Ovarian Cancer'],
      [50,50,'Bladder Cancer'],
      [51,51,'Renal Cancer'],
      [70,70,'Glioblastoma'],
      [106,106,'Gastric Cancer'],
      [130,130,'Malignant Melanoma'],
      [161,161,'Prostate Cancer'],
      [211,211,'Breast Cancer'],
      [229,229,'Lung Cancer'],
      [277,277,'Colon Cancer'],
      [637,637,'Colorectal Cancer'],
      [805,805,'pancreatic adenocarcinoma '],
      [1066,1066,'Hepatocellular Carcinoma'],
    ]
  },
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15},},
  grid: { containLabel: true },
  xAxis: { type: 'category',axisLabel: {fontFamily: "Montserrat",rotate:45,fontSize:8},},
  yAxis: { name: 'Count',axisLabel: {fontFamily: "Montserrat"}},
  grid: {
      left: '5%',
      right: '5%',
      bottom: '12%',
      containLabel: true
  },
  visualMap: {
    orient: 'horizontal',
    left: 'center',
    min: 10,
    max: 100,
    textStyle:{fontFamily: 'Montserrat',fontSize:12},
    text: [' ', 'Count'],
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
          position: 'top',
          valueAnimation: true,
          fontFamily: 'Montserrat'
        },
      encode: {
        // Map the "amount" column to X axis.
        x: 'product',
        // Map the "product" column to Y axis
        y: 'amount'
      }
    }
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>

</body>
</html>