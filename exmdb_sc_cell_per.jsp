<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">






</head>



<body>

<div id="container" style="height: 500px"></div>

<script type="text/javascript" src="echart/echarts.min.js"></script>

<script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};
var option;

<% 

String cancername = "GSE142213";  
if(request.getParameter("cancername")!=null&&!request.getParameter("cancername").equals("")){	
	cancername = request.getParameter("cancername");
}

String path = request.getRealPath("/").replace("\\","/");   //get root real path 
String immpath  = path+"sc_cell_per/"+cancername+"_table.txt";

String cells = "";
int wh_i = 0;
BufferedReader br = new fcon(immpath).getBr();
String str = br.readLine(); // igron 1 line

String title [] = str.split("\t");



StringBuffer sb = new StringBuffer();
StringBuffer echart_in = new StringBuffer();

for(int title_i=0;title_i<title.length;title_i++){
	
	if(title_i == 0 ){
		sb.append(title[title_i]);
		sb.append("=['");
	}else if (title_i == title.length-1 ){
		sb.append(title[title_i]);
	}else{
		sb.append(title[title_i]);
		sb.append("'");
		sb.append(",");
		sb.append("'");
	}
}
sb.append("']"+"\n");

str = br.readLine();

while(str!=null && str.length()>0) {
	wh_i = wh_i + 1;
	String cellarray [] = str.split("\t");
	for(int index_i=0;index_i<cellarray.length;index_i++){
		sb.append(cellarray[index_i]);
		if(index_i == 0 ){
			sb.append("=[");
			echart_in.append(  "{name: '" + cellarray[index_i] + "',type: 'bar',stack: 'total',emphasis: {focus: 'series'},data:"+ cellarray[index_i] + "},\n");
		}else{
			sb.append(",");
		}
		
		
		
	}
	sb.append("]"+"\n");
	str = br.readLine();
}

cells = sb.toString();
out.println(cells);
// System.out.println(cells);
// System.out.println(echart_in.toString());
%>

option = {
              title: {
                  text: 'Data ID : <%=cancername%>',
                  textStyle:{fontFamily: 'Montserrat'},
              },
			  
   toolbox: {
	 feature: {
		magicType: {type: ['stack']},
		saveAsImage: {pixelRatio: 2},
		dataView: {}
		}
  },
  color:["#1B2E5A","#254399","#1A63A0","#21B2D1","#F1EA36","#F4AB20","#D73450","#E14D6B","#B85D9F"],
  tooltip: {
    trigger: 'axis',
    textStyle:{fontFamily: 'Montserrat',fontSize:15},
    axisPointer: {
      // Use axis to trigger tooltip
      type: 'shadow' // 'shadow' as default; can also be 'line' or 'shadow'
    }
  },
  legend: {
	  textStyle:{fontFamily: 'Montserrat',fontSize:15},
	  left: 'center',
	  bottom: 0,
	  type: 'scroll'
  },
  grid: {
	top:'15%',
    left: '0%',
    right: '0%',
    bottom: '10%',
    containLabel: true
  },
  xAxis: {
	max:1,
    type: 'value'
  },
  yAxis: {
      type: 'category',
      data: mix,
      axisLabel: {
          rotate: 45,
          fontFamily: 'Montserrat'
      }
  },
  series: [
			<%=echart_in.toString()%>
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>

 

</body>
</html>