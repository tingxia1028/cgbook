<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/gitalk@1/dist/gitalk.css">
<script src="https://cdn.jsdelivr.net/npm/gitalk@1/dist/gitalk.min.js"></script>
<div id="gitalk-container"></div>
<script>
  var gitalk = new Gitalk({
    "clientID": "59e13207ca973ad27fc7",
    "clientSecret": "71672536db37447deccc73cb06df6d43fa29142a",
    "repo": "cgbook",
    "owner": "TingXia",
    "admin": ["TingXia"],
    "id": location.pathname,      
    "distractionFreeMode": false  
  });
  gitalk.render("gitalk-container");
</script>