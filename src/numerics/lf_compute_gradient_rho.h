    RealType tmpgradx_rho = 0;
    RealType tmpgrady_rho = 0;
    sum_alphax = 0;
    sum_alphay = 0;
    // m_rho[i]=0;
    rho_grady_top  [cell_ooo] = 0;
    rho_grady_bot  [cell_ooo] = 0;
    rho_gradx_right[cell_ooo] = 0;
    rho_gradx_left [cell_ooo] = 0;

    for(int j = 0; j < nb_mat; j++){
      //
      RealType in_rho_k_ooo  = in_rho_k[j][cell_ooo];  
      RealType in_rho_k_m1m1 = in_rho_k[j][cell_m1m1];
      RealType in_rho_k_p1p1 = in_rho_k[j][cell_p1p1];
      RealType in_rho_k_m1p1 = in_rho_k[j][cell_m1p1];
      RealType in_rho_k_p1m1 = in_rho_k[j][cell_p1m1];
      RealType in_rho_k_om1  = in_rho_k[j][cell_om1];
      RealType in_rho_k_op1  = in_rho_k[j][cell_op1];
      RealType in_rho_k_m1o  = in_rho_k[j][cell_m1o];
      RealType in_rho_k_p1o  = in_rho_k[j][cell_p1o];

      zpp = 1.;
      zpm = 1.;
      zmp = 1.;
      zmm = 1.;

      gradx = twelth * (in_rho_k_p1p1 - in_rho_k_m1p1) +
	      third  * (in_rho_k_p1o  -in_rho_k_m1o  ) +
	      twelth * (in_rho_k_p1m1 - in_rho_k_m1m1);

      grady = twelth * (in_rho_k_p1p1 - in_rho_k_p1m1) +
	      third  * (in_rho_k_op1  - in_rho_k_om1 ) +
	      twelth * (in_rho_k_m1p1 - in_rho_k_m1m1);


      grad_tmp = 0.5 * grady + 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zpp     = std::max(in_rho_k_ooo,in_rho_k_p1o);
	  zpp_tmp = std::max(in_rho_k_op1,in_rho_k_p1p1);
	  zpp     = std::max(zpp,zpp_tmp);
	  //zsup=1.;
	  zpp = (zpp-in_rho_k_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpp     = std::min(in_rho_k_ooo,in_rho_k_p1o);
	  zpp_tmp = std::min(in_rho_k_op1,in_rho_k_p1p1);
	  zpp     = std::min(zpp,zpp_tmp);
	  //zsup=0;
	  zpp = (zpp-in_rho_k_ooo) / (grad_tmp);
	}
      grad_tmp = -0.5 * grady + 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zpm     = std::max(in_rho_k_ooo,in_rho_k_p1o);
	  zpm_tmp = std::max(in_rho_k_om1,in_rho_k_p1m1);
	  zpm     = std::max(zpm,zpm_tmp);
	  //zsup=1.;
	  zpm = (zpm-in_rho_k_ooo)/(grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpm      = std::min(in_rho_k_ooo,in_rho_k_p1o);
	  zpm_tmp  = std::min(in_rho_k_om1,in_rho_k_p1m1);
	  zpm      = std::min(zpm,zpm_tmp);
	  //zsup=0;
	  zpm = (zpm - in_rho_k_ooo) / (grad_tmp);
	}

      grad_tmp = 0.5 * grady - 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zmp     = std::max(in_rho_k_ooo,in_rho_k_m1o);
	  zmp_tmp = std::max(in_rho_k_op1,in_rho_k_m1p1);
	  zmp     = std::max(zmp,zmp_tmp);
	  //zsup=1.;
	  zmp = (zmp - in_rho_k_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmp     = std::min(in_rho_k_ooo,in_rho_k_m1o);
	  zmp_tmp = std::min(in_rho_k_op1,in_rho_k_m1p1);
	  zmp     = std::min(zmp,zmp_tmp);
	  //zsup=0;
	  zmp = (zmp - in_rho_k_ooo) / (grad_tmp);
	}

      grad_tmp = -0.5 * grady - 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zmm     = std::max(in_rho_k_ooo,in_rho_k_m1o);
	  zmm_tmp = std::max(in_rho_k_om1,in_rho_k_m1m1);
	  zmm     = std::max(zmm,zmm_tmp);
	  //zsup=1.;
	  zmm = (zmm - in_rho_k_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmm     = std::min(in_rho_k_ooo,in_rho_k_m1o);
	  zmm_tmp = std::min(in_rho_k_om1,in_rho_k_m1m1);
	  zmm     = std::min(zmm,zmm_tmp);
	  //zsup=0;
	  zmm = (zmm - in_rho_k_ooo) / (grad_tmp);
	}

      tmp = std::min(zpp,zmm);
      phi = std::min(zmp,zpm);
      phi = std::min(phi,tmp);
      phi = std::min(0.5*phi,1.0);
      rhok_grady[j][cell_ooo] = phi * grady * h_y;
      rhok_gradx[j][cell_ooo] = phi * gradx * h_x;
      rho_grady_top[cell_ooo] = rho_grady_top[cell_ooo] + rhok_grady[j][cell_ooo]
	* in_c_k[j][cell_ooo]
	+ alphak_grady_top[j][cell_ooo]
	* in_rho_k_ooo
	+ 0.5 * dy * alphak_grady_top[j][cell_ooo]
	* rhok_grady[j][cell_ooo];
      
      rho_grady_bot[cell_ooo] = rho_grady_bot[cell_ooo] + rhok_grady[j][cell_ooo]
	* in_c_k[j][cell_ooo]
	+ alphak_grady_bot[j][cell_ooo]
	* in_rho_k_ooo
	- 0.5 * dy * alphak_grady_bot[j][cell_ooo]
	* rhok_grady[j][cell_ooo];
      
      rho_gradx_right[cell_ooo] = rho_gradx_right[cell_ooo] + rhok_gradx[j][cell_ooo]
	* in_c_k[j][cell_ooo]
	+ alphak_gradx_right[j][cell_ooo]
	* in_rho_k_ooo
	+ 0.5 * dx * alphak_gradx_right[j][cell_ooo]
	*rhok_gradx[j][cell_ooo];

      rho_gradx_left[cell_ooo]= rho_gradx_left[cell_ooo] + rhok_gradx[j][cell_ooo]
	* in_c_k[j][cell_ooo]
	+ alphak_gradx_left[j][cell_ooo]
	* in_rho_k_ooo
	-0.5 * dx * alphak_gradx_left[j][cell_ooo]
	* rhok_gradx[j][cell_ooo];
    }
