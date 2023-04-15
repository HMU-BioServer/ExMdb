<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!DOCTYPE html>
<html style="height: 100%">
<head>
<meta charset="utf-8">

<link rel="stylesheet" href="assets/css/style.css">

<!-- Datatable -->
<link rel="stylesheet" href="assets/css/style.css">
<link rel="stylesheet" type="text/css" href="Datatable/css/table.css">
<link rel="stylesheet" type="text/css" href="Datatable/css/buttons.dataTables.min2.css">
<link rel="stylesheet" href="Datatable/css/reset.css">
<link rel="stylesheet" href="Datatable/css/icon.css">
</head>
<body style="height: 100%; margin:0">

<div id="container" style="width:48%;min-width:500px;height: 400px;float:left"></div>


<script type="text/javascript" src="echart/echarts.min.js"></script>
<script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};
var option;

<%
String GSEID = "GSE110271";
int i = 0;
String path = request.getRealPath("/").replace("\\","/");   //get root real path 

if(request.getParameter("GSEID")!=null&&!request.getParameter("GSEID").equals("")){
	GSEID = request.getParameter("GSEID");
}

String searchname = "";
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}

String items = "200";
if(request.getParameter("items")!=null&&!request.getParameter("items").equals("")){	
	items = request.getParameter("items");
}

String Sql_volcano = "select * from limma_res where i = \"" + GSEID + "\"";
String Sql_volcano_table = "select * from limma_res where i = \"" + GSEID + "\" limit "+ items;

String a = dbhello.ExMdb_Volcano(Sql_volcano);
// System.out.println("finish");
String limma_table = dbhello.limma_table(Sql_volcano_table);
// System.out.println("finish 1");
String barset = dbhello.limma_bar(Sql_volcano);
// System.out.println("finish 2");
// String top_mirna = dbhello.top_mirna(Sql_volcano);
// System.out.println("finish 3");
// String limma_table = "";
// String barset = "";
// String top_mirna = "";
// System.out.println(a);
// System.out.println(BarSet);
// out.println("Nosig="+Nosig);
// out.println("Up="+Up);
// out.println("Down="+Down);



String locxy = "";
locxy = dbhello.exmdb_getxy(searchname,GSEID);


out.println(a);
%>



option = {
  title: {text: 'Volcano plot',
	 	  subtext: '-Log10(p) (Inf)= 400',
		   textStyle: {
			   fontSize: 15,
			   fontFamily: 'Montserrat'
		   }},
  xAxis: {
	nameTextStyle: {fontFamily: "Montserrat"},
    axisLabel: {fontFamily: "Montserrat"},
    scale: true,
    name:"Log2(FC)",
    nameLocation:"middle",//middle start end
    nameGap:30,
  },
  yAxis: {
	nameTextStyle: {fontFamily: "Montserrat",padding: 20},
	axisLabel: {fontFamily: "Montserrat"},
    scale: true,
    name:"-Log10(P)",
    nameLocation:"middle",//middle start end
  },
  toolbox: {
    feature: {
      dataView: { show: true, readOnly: false },
      restore: { show: true },
      saveAsImage: { show: true }
    }
  },
      tooltip: {
    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        trigger: 'item',
        axisPointer: {
          type: 'cross'
        },
        formatter: function (params) { 
          return params.marker+"Gene Name : "+params.data[2]+"<br/>"+params.marker+"Log2(FC) : "+params.data[0] +"<br/>"+params.marker+"-Log10(p-Value) : "+params.data[1];
        }
      },
    color:['#737373','red','green',"#b3ccff"],
    legend: {
    textStyle:{fontFamily: 'Montserrat',fontSize:15},
    left: 'center'
  },
  series: [

    {
      name: 'No Sig',
      type: 'scatter',
      itemStyle: {
      normal: {
       borderWidth: 0.5,
       borderColor: '#fff',
       color: "#cccccc"
      }},
      // prettier-ignore
      data: Nosig
    },
    {
      name: 'Up',
      type: 'scatter',
      emphasis: {focus: 'self'},
      itemStyle: {
      normal: {
       borderWidth: 0.25,
       borderColor: '#fff',
       color: "#F5D043"
                }},
      // prettier-ignore
      data: Up
      },
    {
      name: 'Down',
      type: 'scatter',
      emphasis: {focus: 'self'},
      itemStyle: {
      normal: {
       borderWidth: 0.25,
       borderColor: '#fff',
       color: "#114182"
                }},
      // prettier-ignore
      data: Down    
     },
	  {
         type: 'effectScatter',
         symbolSize: 20,
         color: "#ff6666",
         data: [
			<%=locxy%>
         ]
       },
 	  {
         type: 'scatter',
         emphasis: {focus: 'self'},
         color: "white",
         itemStyle: {
             normal: {
              borderWidth: 0.25,
              borderColor: '#fff',
              color: "white"
                       }},
           data: [
  			<%=locxy%>
           ]
         },
  ]
};



if (option && typeof option === 'object') {
    myChart.setOption(option);
};

myChart.on('click', function (params) {
	if(params.data[2]){
		window.open("exmdb_volcano_out.jsp?searchname="+params.data[2]+"&GSEID="+"<%=GSEID%>")
	}
	
})

</script>
<div id="Bar" style="height: 400px;width:50%;float:left"></div>
<script type="text/javascript">
var dom = document.getElementById("Bar");
var myChart = echarts.init(dom);
var app = {};

var option;



option = {
title: {text: 'Top 10 most significant genes',
	    textStyle: {
		   fontSize: 15,
		   fontFamily: 'Montserrat'
	   },
	},
  legend: {textStyle:{fontFamily: 'Montserrat',fontSize:15}},
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15}},
  dataset: {
	  <%out.println(barset);%>
  },
  xAxis: {
	nameTextStyle: {fontFamily: "Montserrat"},
    axisLabel: {
        fontFamily: 'Montserrat',
        rotate: 30
      },
    type: 'category',
  },
  yAxis: {
	nameTextStyle: {fontFamily: "Montserrat"},
	axisLabel: {fontFamily: "Montserrat"},
  },
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

myChart.on('click', function (params) {

	if(params.name){
		window.open("exmdb_volcano_out.jsp?searchname="+params.name+"&GSEID="+"<%=GSEID%>")
	}
})


</script>
        


<div id="content-d" class="content light" style="padding:2rem;width:100%;height:800px" >


		<div class="row scroll-pane">
		  <div id="admin" class="col s12">
		    <div class="card material-table">
		      <div class="table-header" id="buttons">
		      
		        <div class="actions">
		          <a class="buttons-toggle waves-effect btn-flat nopadding"   href="javascript:;" οnclick="redirect($(this))" val="#"><i class="material-icons">folder_shared</i></a>
		          <a class="search-toggle waves-effect btn-flat nopadding" href="javascript:;" οnclick="redirect($(this))" val="#"><i class="material-icons">search</i></a>
		        </div>
		        
		        <form action="exmdb_volcano.jsp"  method="post" role="form"  style="display: contents;">
				<label style="font-size: 20px;">Items: </label>
				<select name="items" id="items" style="display: block;width: 20%;margin:0 10px;">
				  <option value="200">200</option>
				  <option value="500">500</option>
				  <option value="1000">1000</option>
				  <option value="2000">2000</option>
				</select>
				
				<input style="display:none" type="text"  name="searchname" id="searchname" value="<%=searchname%>">
				<input style="display:none" type="text"  name="GSEID" id="GSEID"  value="<%=GSEID%>">
				
				<button type="submit" class="dt-button buttons-excel buttons-html5">Submit</button>
				</form>
		      </div>
		      <table id="datatable1">
		        <thead>
		         	 <tr>
		         	 			<th rowspan="1" colspan="1" style="width:60px;">Gene Name</th> 
		         	 			<th rowspan="1" colspan="1" style="width:60px;">Control Mean</th>
								<th rowspan="1" colspan="1" style="width:60px;">Case Mean</th>
								<th rowspan="1" colspan="1" style="width:60px;">-Log10(p)</th>	
								<th rowspan="1" colspan="1" style="width:60px;">Log2(FC)</th>
		            </tr> 
		        </thead>
		        <tbody>
		        <%out.println(limma_table);%>
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

</body>
</html>
    