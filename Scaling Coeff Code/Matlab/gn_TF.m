function TF = gn_TF(p_lead,p_TF,polate_mode)

tip2can_T.mag = [403	212	163	156	157	159	159	159	150	154	154	149	147	139	134	127	119	110	100	92	83	73	65	58	51	46	44	49	55	61	71	62	88	96	108	118	128	132	140];
tip2can_T.phase = [-83	-106	-121	-134	-148	-152	-159	-164	-165	-170	-173	-176	-177	178	175	172	168	166	161	157	153	146	139	128	116	98	81	61	53	42	32	27	23	20	16	14	14	12	13];
tip2can_T.x = (0:1:length(tip2can_T.mag)-1);
% figure
% subplot(2,1,1)
% plot(tip2can_T.x,tip2can_T.mag)
% subplot(2,1,2)
% plot(tip2can_T.x,tip2can_T.phase)

tip2can_V.mag = [1.12	1.25	1.38	1.475	1.515	1.552	1.554	1.57	1.575	1.55	1.52	1.5	1.45	1.42	1.37	1.31	1.26	1.225	1.19	1.18	1.16	1.13	1.12	1.1	1.08	1.06	1.03	1.03	1.03	1.03	1.04	1.08	1.16	1.24	1.3	1.4	1.49	1.57	1.66	1.75];
tip2can_V.phase = [-9	-3	1	3.8	6	8	9	11	13	16	17.6	20	23	26	29	33	38	42	43	46	48	52	54	56	60	66	69	74	78	87	96	107	118	122	129	136	142	146	150	155];
tip2can_V.x = (0:1:length(tip2can_V.mag)-1);
% figure
% subplot(2,1,1)
% plot(tip2can_V.x,tip2can_V.mag)
% subplot(2,1,2)
% plot(tip2can_V.x,tip2can_V.phase)

ring2can_T.mag = [128	96	91	89	87	88	88	87	87	86	83	82	77	76	71	65	60	54	52	44	42	34	27	25	22	25	23	27	34	38	43	49	54	59	64	69	76	77];
ring2can_T.phase = [-112	-122	-144	-148	-150	-153	-155	-160	-164	-166	-169	-175	-176	-178	-179	177	171	168	164	160	150	140	138	119	99	85	74	51	48	39	33	29	25	22	25	18	21	18];
ring2can_T.x = (0:1:length(ring2can_T.mag)-1);
% figure
% subplot(2,1,1)
% plot(ring2can_T.x,ring2can_T.mag)
% subplot(2,1,2)
% plot(ring2can_T.x,ring2can_T.phase)


ring2can_V.mag = [1.04	1.18	1.28	1.35	1.42	1.45	1.48	1.48	1.5	1.52	1.52	1.53	1.53	1.52	1.51	1.49	1.49	1.46	1.42	1.37	1.33	1.28	1.25	1.18	1.18	1.15	1.1	1.06	1.03	1.01	1.01	1.03	1.05	1.12	1.16	1.2	1.3	1.4	1.51	1.59];
ring2can_V.phase = [-10	-5	-2	1	3	4	6	7	8	10	10	12	13	14	16	18	20	22	25	28	31	36	38	42	45	47	53	60	69	79	88	98	104	115	120	126	132	139	147	149];
ring2can_V.x = (0:1:length(ring2can_V.mag)-1);
% figure
% subplot(2,1,1)
% plot(ring2can_V.x,ring2can_V.mag)
% subplot(2,1,2)
% plot(ring2can_V.x,ring2can_V.phase)

switch p_lead
    case 'ring'
        switch p_TF
            case 'T'
                TF = ring2can_T;
            case 'V'
                TF = ring2can_V;
            otherwise
                TF = [];
                disp('The optional parameter is ''T'' and ''V''.');
                return
        end
        
    case 'tip'
        switch p_TF
            case 'T'
                TF = tip2can_T;
            case 'V'
                TF = tip2can_V;
            otherwise
                TF = [];
                disp('The optional parameter is ''T'' and ''V''.');
                return
        end
    otherwise
        TF = [];
        disp('The optional parameter is ''ring'' and ''tip''.');
        return
end

switch polate_mode
    case 'interp'
        label.x = TF.x(end);
        label.y1 = TF.mag(end);
        label.y2 = TF.phase(end);
        TF.x = [TF.x 43];
        TF.mag = [TF.mag 0]; TF.phase = [TF.phase TF.phase(end)];
        TF.mag = interp1(TF.x,TF.mag,(0:1:43));
        TF.phase = interp1(TF.x,TF.phase,(0:1:43));
        TF.x = (0:1:43);
        figure
        subplot(2,1,1)
        plot(TF.x,TF.mag)
        hold on
        scatter(label.x,label.y1);
        subplot(2,1,2)
        plot(TF.x,TF.phase);
        hold on
        scatter(label.x,label.y2);
        disp('The interpolation function is used.');
    case 'extrap'
        label.x = TF.x(end);
        label.y1 = TF.mag(end);
        label.y2 = TF.phase(end);
        Test.mag = TF.mag;
        Test.phase=TF.phase;
        TF.mag = interp1(TF.x,TF.mag,(0:1:43),'pchip','extrap');
        TF.phase = interp1(TF.x,TF.phase,(0:1:43),'pchip','extrap');
        TF.re = interp1(TF.x,real(Test.mag.*exp(1j*Test.phase./180*pi)),(0:1:43),'pchip','extrap');
        TF.im = interp1(TF.x,imag(Test.mag.*exp(1j*Test.phase./180*pi)),(0:1:43),'pchip','extrap');
        TF.x = (0:1:43);
        figure
        subplot(2,1,1)
        plot(TF.x,TF.mag,'b')
        title(['Magnitude of ', p_TF,' transfer function of ',p_lead])
        hold on
        plot(TF.x,abs(TF.re+1j*TF.im),'r')
        scatter(label.x,label.y1,'k');
        subplot(2,1,2)
        plot(TF.x,TF.phase,'b');
        hold on
        plot(TF.x,angle(TF.re+1j*TF.im)./pi*180,'r')
        title(['Phase of ', p_TF,' transfer function of ',p_lead])
        scatter(label.x,label.y2,'k');
        disp('The exterpolation function is used.');
        
    case 'no'
        figure
        subplot(2,1,1)
        plot(ring2can_V.x,ring2can_V.mag)
        subplot(2,1,2)
        plot(ring2can_V.x,ring2can_V.phase)
        disp('The interpolation/exterpolation is not performed, raw TF is used.');
    otherwise 
        disp('The parameter is not available.');
        return;
end

TF.y = TF.mag.*exp(1j*TF.phase*pi/180);

           










