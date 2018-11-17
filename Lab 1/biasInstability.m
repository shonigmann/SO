function BI = biasInstability(t,W,Tbi)

    X = t;
    X(1) = W(1);
    
    for k=1:length(t)-1

        if(mod(t(k+1),Tbi) == 0)
            X(k+1) = W(k+1);
        else
            X(k+1) = X(k);
        end
    end
    
    figure(11);
    plot(t,X);
    title('Generated Bias Instability Signal Created from Unity Variance White Noise');
    xlabel('Time, s');
    ylabel('Signal Magnitude');
    
    BI=X;
end