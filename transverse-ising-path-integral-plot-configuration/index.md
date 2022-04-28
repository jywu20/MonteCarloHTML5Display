在这个页面上，你可以设置二维正方晶格上的伊辛模型
$$
H = - J \sum_{\lang \boldsymbol{i}, \boldsymbol{j} \rang} \sigma_{\boldsymbol{i}} \sigma_{\boldsymbol{j}} - h \sum_{\boldsymbol{i}} \sigma_{\boldsymbol{i}}
$$
的磁场$h$（我们固定$J=1$）和系统温度$T$，以及系统大小$L$，观察自旋构型会发生什么变化。提示：$J = 1$时，$T_\text{c} = 2.2691853142$。观察接近临界点和远离临界点时的系统行为。

<canvas id="myCanvas" width="600" height="600"></canvas>
<form>
  <label for="magnetic_field">$h$</label>
  <input type="text" id="magnetic_field" name="magnetic_field"><br>
  <label for="temperature">$T$</label>
  <input type="text" id="temperature" name="temperature"><br>
  <label for="lattice_size">$L$</label>
  <input type="text" id="lattice_size" name="lattice_size"><br>
  <button id="start_show_config">开始运行模拟</button>
</form>

<script src="math-js-10.4.0.js"></script>
<script src="utils.js"></script>
<script src="ising-core.js"></script>
<script src="lattice.js"></script>
<script src="app.js"></script>