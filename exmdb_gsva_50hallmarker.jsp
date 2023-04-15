<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">






</head>



<body>
<img src="gsva_50hallmarker/img/gsva_left.png" style="position: absolute;width:28.5%;margin:5% 0 0 0">
<div id="container" style="height: 750px;"></div>

<script type="text/javascript" src="echart/echarts.min.js"></script>

<script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};
var option;


<% 

String cancername = "GSE125442";  
if(request.getParameter("cancername")!=null&&!request.getParameter("cancername").equals("")){	
	cancername = request.getParameter("cancername");
}

String path = request.getRealPath("/").replace("\\","/");   //get root real path 
String immpath  = path+"gsva_50hallmarker/"+cancername+"_GSVA.txt";

String name = "";
String data = "";
int wh_i = 0;
BufferedReader br = new fcon(immpath).getBr();
String str = br.readLine();

String title [] = str.split("\t");



StringBuffer sb_1 = new StringBuffer();
StringBuffer sb_2 = new StringBuffer();

while(str!=null && str.length()>0) {
	String cellarray [] = str.split("\t");
	sb_1.append("'"+cellarray[0]+"',");
	if(Double.valueOf(cellarray[1]) > 0){sb_2.append("{value:"+cellarray[1]+",label: labelLeft},");}
	else{sb_2.append("{value:"+cellarray[1]+",label: labelRight},");}
	
	str = br.readLine();
}

name = sb_1.toString();
data = sb_2.toString();
// out.println(name);
// out.println(data);

// System.out.println(name);
// System.out.println(data);

%>




const labelRight = {
	position: 'right'
};

const labelLeft = {
	position: 'left'
};
		
option = {
  title: {
    text: 'GSVA: 50 Hallmarker',
    subtext: 'Cancer vs Normal',
    textStyle:{fontFamily: 'Montserrat',fontSize:30},
  },
  tooltip: {
    trigger: 'axis',
    textStyle:{fontFamily: 'Montserrat',fontSize:15},
    axisPointer: {
      type: 'shadow'
    }
  },
  grid: {
		top:'6%',
	    left: '30%',
	    right: '4%',
	    bottom: '6%',
	    containLabel: true
	  },
  xAxis: {
    type: 'value',
    position: 'top',
    splitLine: {
      lineStyle: {
        type: 'dashed'
      }
    }
  },
  yAxis: {
    type: 'category',
    axisLine: { show: false },
    axisLabel: { show: false },
    axisTick: { show: false },
    splitLine: { show: false },
    data: [
    	<%=name%>
    ]
  },
  visualMap: {
    orient: 'horizontal',
    left: 'center',
    min: -3,
    max: 3,
    text: ['Cancer', 'Normal'],
    textStyle:{fontFamily: 'Montserrat',fontSize:12},
    // Map the score column to color
    dimension: 0,
    inRange: {
    	 color: ['#F5D043','#114182']
    }
  },
  series: [
    {
      name: 't score',
      type: 'bar',
      stack: 'Total',
      label: {
        show: true,
        formatter: '{b}',
        textStyle:{fontFamily: 'Montserrat',fontSize:12},
        
      },
      data: [
    	  <%=data%>
      ]
    }
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

        </script>

 

</body>
</html>