
var a=0
$(document).ready(function(){
    $(".sidebar-toggle").click(function(){

        $("span.logo").toggle(speed="slow",ease="linear");
         /* $("this").toggle("slide", { direction: "right" }, 1000px);*/
         if(a==0){
       $(".sidebar-toggle").css({"margin-left":"-232px"});
       $(".col-sm-3").hide();
          $(".col-sm-8").css({"width":"90%"});
          $(".row").css({"width":"100%"})
       a=1
   }
   else{
   	 $(".sidebar-toggle").css({"margin-left":"-10px"});
   	 $(".col-sm-3").show();
   	
          $(".col-sm-8").css({"width":" 66.66666667%"});
           $(".row").css({"width":"80%"})
   	
   	a=0
   }
         
       
    });
});
