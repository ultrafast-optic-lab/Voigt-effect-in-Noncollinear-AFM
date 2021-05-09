clc;
clear;

s_d = 3.015 + 16.75 * 1i;
d_group = [0 -0.0005 * s_d];
Gxxxx_group = [0.004 * s_d 0.00605 * s_d];
Gxxzz_group = [0.002 * s_d 0.00205 * s_d];
Gyyxx_group = [0.003 * s_d 0.00301 * s_d];
Gxxyy = 0.00405 * s_d;
Gzyzy = 0.00100 * s_d;
Gyyyy = 0.00602 * s_d;
reslut = [];

for m = 1:2
    for n = 1:2
        d = d_group(m);
        Gxxxx = Gxxxx_group(n);
        Gxxzz = Gxxzz_group(n);
        Gyyxx = Gyyxx_group(n);

        Mx1 = sqrt(3)/2; My1 = 0; Mz1 = 1/2;
        sxx1 = s_d/3 + Gxxxx * (Mx1)^2 + Gxxyy * (My1)^2 + Gxxzz*(Mz1)^2;
        syy1 = s_d/3 + d/3 + Gyyxx * (Mx1)^2 + Gyyyy * (My1)^2 + Gyyxx * (Mz1)^2;
        szz1 = s_d/3 + Gxxzz * (Mx1)^2 + Gxxyy * (My1)^2 + Gxxxx * (Mz1)^2;
        sxy1 = 2 * Gzyzy * Mx1 * My1; syx1 = 2 * Gzyzy * Mx1 * My1;
        sxz1 = (Gxxxx - Gxxyy) * Mx1 * Mz1; szx1 = (Gxxxx - Gxxyy) * Mx1 * Mz1;
        szy1 = 2 * Gzyzy * My1 * Mz1; syz1 = 2 * Gzyzy * My1 * Mz1;
        s1 = [sxx1    sxy1    sxz1
              syx1    syy1    syz1
              szx1    szy1    szz1];

        Mx2 = 0; My2 = 0; Mz2=-1;
        sxx2 = s_d/3 + Gxxxx * (Mx2)^2 + Gxxyy * (My2)^2 + Gxxzz * (Mz2)^2;
        syy2 = s_d/3 + d/3 + Gyyxx * (Mx2)^2 + Gyyyy * (My2)^2 + Gyyxx * (Mz2)^2;
        szz2 = s_d/3 + Gxxzz * (Mx2)^2 + Gxxyy * (My2)^2 + Gxxxx * (Mz2)^2;
        sxy2 = 2 * Gzyzy * Mx2 * My2; syx2 = 2 * Gzyzy * Mx2 * My2;
        sxz2 = (Gxxxx - Gxxyy) * Mx2 * Mz2; szx2 = (Gxxxx - Gxxyy) * Mx2 * Mz2;
        szy2 = 2 * Gzyzy * My2 * Mz2; syz2 = 2 * Gzyzy * My2 * Mz2;
        s2 = [sxx2    sxy2    sxz2
              syx2    syy2    syz2
              szx2   szy2    szz2];

        Mx3 = -sqrt(3)/2; My3 = 0; Mz3 = 1/2;
        sxx3 = s_d/3 + Gxxxx * (Mx3)^2 + Gxxyy * (My3)^2 + Gxxzz * (Mz3)^2;
        syy3 = s_d/3 + d/3 + Gyyxx * (Mx3)^2 + Gyyyy * (My3)^2 + Gyyxx*(Mz3)^2;
        szz3 = s_d/3 + Gxxzz * (Mx3)^2 + Gxxyy * (My3)^2 + Gxxxx * (Mz3)^2;
        sxy3 = 2 * Gzyzy * Mx3 * My3; syx3 = 2 * Gzyzy * Mx3 * My3;
        sxz3 = (Gxxxx - Gxxyy) * Mx3 * Mz3; szx3 = (Gxxxx - Gxxyy) * Mx3 * Mz3;
        szy3 = 2 * Gzyzy * My3 * Mz3; syz3 = 2 * Gzyzy * My3 * Mz3;
        s3 = [sxx3    sxy3    sxz3
              syx3    syy3    syz3
              szx3   szy3    szz3];     
        s = s1 + s2 + s3;

        polarization_rotation = [];
        qq1 = []; qq3 = []; qq2 = []; qq4 = [];
        qqe1 = []; qqe3 = []; qqe2 = []; qqe4 = [];
        for u = [1:1:89 91:1:179 181:1:269 271:1:359]
            uu = u * pi / 180;
            R = [cos(uu) sin(uu) 0.0
                -sin(uu) cos(uu) 0.0
                    0.0     0.0  1.0];
            s_R = zeros(3,3);
            
            for i = 1:3
                for j = 1:3
                    for x = 1:3
                        for y = 1:3
                            s_R(i,j) = s_R(i,j) + R(i,x) * R(j,y) * s(x,y);
                        end
                    end
                end
            end    

            s_0 = 1.0; n_0 = sqrt(s_0);
            q = (0 * pi / 180); ax = sin(q); az = cos(q); nx = n_0*ax;
            sxx = s_R(1,1); sxy = s_R(1,2); sxz = s_R(1,3);
            syx = s_R(2,1); syy = s_R(2,2); syz = s_R(2,3);
            szx = s_R(3,1); szy = s_R(3,2); szz = s_R(3,3);
            syms x
            l1 = 1; l2 = 2; l3 = 3; l4 = 4;
            nz = solve(szz * x^4 + (szx + sxz) * nx * x^3 - [(syy * szz - syz * (szy)) + (szz * sxx - szx * (sxz)) - (sxx + szz) * nx^2] * x^2 + ...
                ((szx + sxz) * (syy + nx^2) - (syz * sxy + szy * syx)) * nx * x + sxx * nx^4 - (sxx * szz - szx * (sxz)) * nx^2 - (sxx * syy -syx * sxy) * nx^2 + ...
                sxx * syy * szz + sxy * syz * szx + syx * szy * sxz - sxx * syz * szy - syy * szx * sxz - szz * sxy * syx);
            if real(nz(l1)) > 0 && real(nz(l2)) > 0
              if imag(nz(l1)) < imag(nz(l2))
                  qq1 = [qq1 nz(l1)]; nz1 = nz(l1);
                  qq3 = [qq3 nz(l2)]; nz3 = nz(l2);
              else
                  qq1 = [qq1 nz(l2)]; nz1 = nz(l2);
                  qq3 = [qq3 nz(l1)]; nz3 = nz(l1);
              end
            elseif real(nz(l1)) > 0 && real(nz(l3)) > 0
              if imag(nz(l1)) < imag(nz(l3))
                  qq1 = [qq1 nz(l1)]; nz1 = nz(l1);
                  qq3 = [qq3 nz(l3)]; nz3 = nz(l3);
              else
                  qq1 = [qq1 nz(l3)]; nz1 = nz(l3);
                  qq3 = [qq3 nz(l1)]; nz3 = nz(l1);
              end
            elseif real(nz(l1)) > 0 && real(nz(l4)) > 0
              if imag(nz(l1)) < imag(nz(l4))
                  qq1 = [qq1 nz(l1)]; nz1 = nz(l1);
                  qq3 = [qq3 nz(l4)]; nz3 = nz(l4);
              else
                  qq1 = [qq1 nz(l4)]; nz1 = nz(l4);
                  qq3 = [qq3 nz(l1)]; nz3 = nz(l1);
              end
            elseif real(nz(l2)) > 0 && real(nz(l3)) > 0
              if imag(nz(l2)) < imag(nz(l3))
                  qq1 = [qq1 nz(l2)]; nz1 = nz(l2);
                  qq3 = [qq3 nz(l3)]; nz3 = nz(l3);
              else
                  qq1 = [qq1 nz(l3)]; nz1 = nz(l3);
                  qq3 = [qq3 nz(l2)]; nz3 = nz(l2);
              end
            elseif real(nz(l2)) > 0 && real(nz(l4)) > 0
              if imag(nz(l2)) < imag(nz(l4))
                  qq1 = [qq1 nz(l2)]; nz1 = nz(l2);
                  qq3 = [qq3 nz(l4)]; nz3 = nz(l4);
              else
                  qq1 = [qq1 nz(l4)]; nz1 = nz(l4);
                  qq3 = [qq3 nz(l2)]; nz3 = nz(l2);
              end
            elseif real(nz(l3)) > 0 && real(nz(l4)) > 0
              if imag(nz(l3)) < imag(nz(l4))
                  qq1 = [qq1 nz(l3)]; nz1 = nz(l3);
                  qq3 = [qq3 nz(l4)]; nz3 = nz(l4);
              else
                  qq1 = [qq1 nz(l4)]; nz1 = nz(l4);
                  qq3 = [qq3 nz(l3)]; nz3 = nz(l3);
              end
            end
            nz2 = 0; nz4 = 0;

            p1 = ((nx^2 - szz) * ((nz1)^2 - sxx) - (nx * nz1 + sxz) * (nx * nz1 + szx)) / ((nx^2 - szz) * sxy + (nx * nz1 + sxz) * szy);
            p2 = ((nx^2 - szz) * ((nz2)^2 - sxx) - (nx * nz2 + sxz) * (nx * nz2 + szx)) / ((nx^2 - szz) * sxy + (nx * nz2 + sxz) * szy);
            p3 = ((nx^2 - szz) * ((nz3)^2 - sxx) - (nx * nz3 + sxz) * (nx * nz3 + szx)) / ((nx^2 - szz) * sxy + (nx * nz3 + sxz) * szy);
            p4 = ((nx^2 - szz) * ((nz4)^2 - sxx) - (nx * nz4 + sxz) * (nx * nz4 + szx)) / ((nx^2 - szz) * sxy + (nx * nz4 + sxz) * szy);
            q1 = (szy * ((nz1)^2 - sxx) + sxy * (nz1 * nx + szx)) / (szy * (nx * nz1 + sxz) + sxy * (nx^2 - szz));
            q2 = (szy * ((nz2)^2 - sxx) + sxy * (nz2 * nx + szx)) / (szy * (nx * nz2 + sxz) + sxy * (nx^2 - szz));
            q3 = (szy * ((nz3)^2 - sxx) + sxy * (nz3 * nx + szx)) / (szy * (nx * nz3 + sxz) + sxy * (nx^2 - szz));
            q4 = (szy * ((nz4)^2 - sxx) + sxy * (nz4 * nx + szx)) / (szy * (nx * nz4 + sxz) + sxy * (nx^2 - szz));
            P1 = nx * (p1 * q3 - p3 * q1) + (nz1 * p3 - nz3 * p1);
            Q1 = p3 - p1;
            S1 = nx * (nz1 * p1 * q3 -nz3 * p3 * q1) + nz1 * nz3 * (p3 - p1);
            T1 = (nz3 * p3 - nz1*p1);
            rank = n_0 * (az)^2 * P1 + (n_0)^2 * az * Q1 + az * S1 + n_0 * T1;
            rss = (n_0 * (az)^2 * P1 + (n_0)^2 * az * Q1 - az * S1 - n_0 * T1)/rank;
            rps = (2 * n_0 * az * (nz1 - nz3 - nx * (q1 - q3))) / rank;
            rsp = (2 * n_0 * az * p1 * p3 * (nz3 - nz1)) / rank;
            rpp = (-n_0 * (az)^2 * P1 + (n_0)^2 * az * Q1 - az * S1 + n_0 * T1) / rank;
            angle_s = -rps / rss;
            angle_p = rsp / rpp;
            polarization_rotation = [polarization_rotation (atan(real(angle_s)) * 10^3)];
        end
        reslut = [reslut ; polarization_rotation];
    end
end

polarization_rotation_deitaG0_d0 = reslut(1,:);
polarization_rotation_deitaG_d0 = reslut(2,:);
polarization_rotation_deitaG0_d = reslut(3,:);
polarization_rotation_deitaG_d = reslut(4,:);

u = [1:1:89 91:1:179 181:1:269 271:1:359];
h = plot(u,real(polarization_rotation_deitaG_d),'b',u,real(polarization_rotation_deitaG_d0),'r',u,real(polarization_rotation_deitaG0_d),'k',u,real(polarization_rotation_deitaG0_d0),'m','LineWidth',4);
title('Simulated polarization rotation')
xlabel('Sample orientation(degree)')
ylabel('polarization rotation(arb. unit)')
legend('deitas0!=0 & deitaG!=0','deitas0=0 & deitaG!=0','deitas0!=0 & deitaG=0','deitas0=0 & deitaG=0')
