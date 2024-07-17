classdef ShaderOperator  < handle
    % Shader constructor Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        imagingMode
        kernel
        kernel1d
        myOperator

        shaderHandleConv1
        shaderHandleConv2
        shaderHandle3
        shaderUniforms
        defaultParams

        texhandle1
        texhandle2
    end
    
    methods

        function initOPENGLShader(this)

            AssertOpenGL;
            InitializeMatlabOpenGL;

            if ismac
                this.imagingMode = kPsychNeedRetinaResolution;
            else
                this.imagingMode = [];
            end

        end

        function checkOPENGLShaderInit(this)

            global GL;
            persistent initialized;

            if isempty(initialized)

                % Make sure GLSL and fragmentshaders are supported on first call:
                AssertGLSL;
            
                % Query supported extensions:
                extensions = glGetString(GL.EXTENSIONS);
                if isempty(findstr(extensions, 'GL_ARB_fragment_shader'))
                    % No fragment shaders: This is a no go!
                    error('Sorry, this function does not work on your graphics hardware due to lack of sufficient support for fragment shaders.');
                end
            
                % Clear any OpenGL error state.
                while (glGetError~=GL.NO_ERROR)
                end
            
                maxuniforms = glGetIntegerv(GL.MAX_FRAGMENT_UNIFORM_COMPONENTS);
                
                % We are initialized:
                initialized = 1;
            end 


        end

        function buildGaussKernel(this,npx,mu,sigma)

            % build 2d kernel, which will be used as a texture
            this.kernel = fspecial("gaussian",npx,sigma);
            this.kernel = (this.kernel./sum(this.kernel,'all'));
            
            % build 1d slice (also normalized to 1), for separable 2d
            % Convolution
            this.kernel1d = normpdf(1:npx,mu,sigma);
            this.kernel1d = (this.kernel1d./sum(this.kernel1d,'all'));

        end

        function replaceGaussKernel(this,npx,mu,newsigma)

            global GL;

            % clamp to reasonable value of SD
            if newsigma == 0
                % impulse response
                this.kernel = zeros(npx,npx);
                this.kernel((npx-1)/2,(npx-1)/2) = 1;
                
                this.kernel1d = zeros(npx,1);
                this.kernel1d((npx-1)/2) = 1;
            else
                this.kernel = fspecial("gaussian",npx,newsigma);
                this.kernel = (this.kernel./sum(this.kernel,'all'));
    
                this.kernel1d = normpdf(1:npx,mu,newsigma);
                this.kernel1d = (this.kernel1d./sum(this.kernel1d,'all'));
            end

            kernelw = length(this.kernel1d);
            kernelh = 1;

            glActiveTexture(GL.TEXTURE1);
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, this.texhandle1);
            glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, GL.LUMINANCE_FLOAT32_APPLE, kernelw, kernelh, 0, GL.LUMINANCE, GL.FLOAT, moglsingle(this.kernel1d));

            % Make sure we use nearest neighbour sampling:
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);

            % And that we clamp to edge:
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP);
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP);

            % Default CLUT setup done: Switch back to texture unit 0:
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
            glActiveTexture(GL.TEXTURE0);

            kernelw = 1;
            kernelh = length(this.kernel1d);

            glActiveTexture(GL.TEXTURE1);
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, this.texhandle2);
            glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, GL.LUMINANCE_FLOAT32_APPLE, kernelw, kernelh, 0, GL.LUMINANCE, GL.FLOAT, moglsingle(this.kernel1d));

            % Make sure we use nearest neighbour sampling:
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);

            % And that we clamp to edge:
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP);
            % glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP);

            % Default CLUT setup done: Switch back to texture unit 0:
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
            glActiveTexture(GL.TEXTURE0);

        end

        function init2dConvGaussShader(this,w)

            global GL;

            kernelw = length(this.kernel1d);
            kernelh = 1;
            hwx = (kernelw - 1) / 2;
            hwy = (kernelh - 1) / 2;

            this.myOperator = CreateGLOperator(w, kPsychNeed32BPCFloat);
            Screen('GetWindowInfo', this.myOperator);

            this.shaderHandleConv1 = LoadGLSLProgramFromFiles(which('Conv2dMatt.frag.txt'), 1);
            
            % Assign proper texture units for input image and clut:
            glUseProgram(this.shaderHandleConv1);
            
            this.shaderUniforms.shaderImage1 = glGetUniformLocation(this.shaderHandleConv1, 'Image');
            glUniform1i(this.shaderUniforms.shaderImage1, 0);
            
            this.shaderUniforms.shaderClut1  = glGetUniformLocation(this.shaderHandleConv1, 'Kernel');
            glUniform1i(this.shaderUniforms.shaderClut1, 1);
            
            this.shaderUniforms.shaderKernelsizeX1  = glGetUniformLocation(this.shaderHandleConv1, 'KernelHalfWidthX');
            glUniform1f(this.shaderUniforms.shaderKernelsizeX1, hwx);
            
            this.shaderUniforms.shaderKernelsizeY1  = glGetUniformLocation(this.shaderHandleConv1, 'KernelHalfWidthY');
            glUniform1f(this.shaderUniforms.shaderKernelsizeY1, hwy);
            
            glUseProgram(0);
            
            glActiveTexture(GL.TEXTURE1);
            this.texhandle1 = glGenTextures(1);
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, this.texhandle1);
            glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, GL.LUMINANCE_FLOAT32_APPLE, kernelw, kernelh, 0, GL.LUMINANCE, GL.FLOAT, moglsingle(this.kernel1d));
            
            % Make sure we use nearest neighbour sampling:
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);
            
            % And that we clamp to edge:
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP);
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP);
            
            % Default CLUT setup done: Switch back to texture unit 0:
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
            glActiveTexture(GL.TEXTURE0);
            configstring = sprintf('TEXTURERECT2D(%i)=%i', 1, this.texhandle1);
            
            count = CountSlotsInGLOperator(this.myOperator);
            if count > 0
                Screen('HookFunction', this.myOperator, 'AppendBuiltin', 'UserDefinedBlit', 'Builtin:FlipFBOs', '');
            end
            
            if count == 0
                % Count was 0, so its now two: Change operator to be at least dual-pass capable:
                Screen('HookFunction', this.myOperator, 'ImagingMode', mor(kPsychNeedDualPass, Screen('HookFunction', this.myOperator, 'ImagingMode')));
            else
                % Change operator to be multi-pass capable:
                Screen('HookFunction', this.myOperator, 'ImagingMode', mor(kPsychNeedMultiPass, Screen('HookFunction', this.myOperator, 'ImagingMode')));
            end
            
            % Add shader to user defined blit chain of the proxy:
            Screen('HookFunction', this.myOperator, 'AppendShader', 'UserDefinedBlit', 'MyShader1', this.shaderHandleConv1, configstring);
            
            % Need a ping-pong op for second convolution pass:
            Screen('HookFunction', this.myOperator, 'AppendBuiltin', 'UserDefinedBlit', 'Builtin:FlipFBOs', '');
            
            
            kernelw = 1;
            kernelh = length(this.kernel1d);
            hwx = (kernelw - 1) / 2;
            hwy = (kernelh - 1) / 2;
            
            this.shaderHandleConv2 = LoadGLSLProgramFromFiles(which('Conv2dMatt.frag.txt'), 1);
            
            % Assign proper texture units for input image and clut:
            glUseProgram(this.shaderHandleConv2);
            
            this.shaderUniforms.shaderImage2 = glGetUniformLocation(this.shaderHandleConv2, 'Image');
            glUniform1i(this.shaderUniforms.shaderImage1, 0);
            
            this.shaderUniforms.shaderClut2  = glGetUniformLocation(this.shaderHandleConv2, 'Kernel');
            glUniform1i(this.shaderUniforms.shaderClut2, 1);
            
            this.shaderUniforms.shaderKernelsizeX2  = glGetUniformLocation(this.shaderHandleConv2, 'KernelHalfWidthX');
            glUniform1f(this.shaderUniforms.shaderKernelsizeX2, hwx);
            
            this.shaderUniforms.shaderKernelsizeY2  = glGetUniformLocation(this.shaderHandleConv2, 'KernelHalfWidthY');
            glUniform1f(this.shaderUniforms.shaderKernelsizeY2, hwy);
            
            glUseProgram(0);
            
            glActiveTexture(GL.TEXTURE1);
            this.texhandle2 = glGenTextures(1);
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, this.texhandle2);
            glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, GL.LUMINANCE_FLOAT32_APPLE, kernelw, kernelh, 0, GL.LUMINANCE, GL.FLOAT, moglsingle(this.kernel1d));
            
            % Make sure we use nearest neighbour sampling:
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);
            
            % And that we clamp to edge:
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP);
            glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP);
            
            % Default CLUT setup done: Switch back to texture unit 0:
            glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
            glActiveTexture(GL.TEXTURE0);
            configstring = sprintf('TEXTURERECT2D(%i)=%i', 1, this.texhandle2);
            
            Screen('HookFunction', this.myOperator, 'AppendShader', 'UserDefinedBlit', 'MyShader2', this.shaderHandleConv2, configstring);
            
            if bitand(Screen('HookFunction', this.myOperator, 'ImagingMode'), mor(kPsychNeed16BPCFloat, kPsychNeed32BPCFloat)) == 0
                % Not yet set. Choose highest precision:
                disp('!')
                Screen('HookFunction', this.myOperator, 'ImagingMode', mor(kPsychNeed32BPCFloat, Screen('HookFunction', this.myOperator, 'ImagingMode')));
                if debug > 3
                    fprintf('Add2DSeparableConvolutionToGLOperator: Increasing precision of operator to 32bpc float.\n');
                end
            end

        end

        function addShader(this, shaderfname, defaults)

            this.defaultParams.gazePosition = [defaults(1),defaults(2)];
            this.defaultParams.radpx = defaults(3);
            this.defaultParams.slope = defaults(4);
            this.defaultParams.nxtoslope = defaults(5);
            this.defaultParams.flatwidth = defaults(6);
            this.defaultParams.coeff = defaults(7); 

            this.defaultParams.rPeriph = defaults(8);
            this.defaultParams.gPeriph = defaults(9);
            this.defaultParams.bPeriph = defaults(10);
            
            this.defaultParams.rFov = defaults(11);
            this.defaultParams.gFov = defaults(12);
            this.defaultParams.bFov = defaults(13);

            this.defaultParams.invgamma = defaults(14);

            Screen('HookFunction', this.myOperator, 'AppendBuiltin', 'UserDefinedBlit', 'Builtin:FlipFBOs', '');

            this.shaderHandle3 = LoadGLSLProgramFromFiles(which(shaderfname), 1);

            % Bind texture unit 0 to shader as input source for the image:
            glUseProgram(this.shaderHandle3);

            this.shaderUniforms.shaderImage3 = glGetUniformLocation(this.shaderHandle3, 'Image');
            glUniform1i(this.shaderUniforms.shaderImage3, 0);

            this.shaderUniforms.gazePosition = glGetUniformLocation(this.shaderHandle3, 'gazePosition');
            glUniform2f(this.shaderUniforms.gazePosition, ...
                this.defaultParams.gazePosition(1), this.defaultParams.gazePosition(2));

            this.shaderUniforms.gazeRadius = glGetUniformLocation(this.shaderHandle3, 'gazeRadius');
            glUniform1f(this.shaderUniforms.gazeRadius, this.defaultParams.radpx);

            this.shaderUniforms.colorScalePeriph = glGetUniformLocation(this.shaderHandle3, 'colorScalePeriph');
            glUniform3f(this.shaderUniforms.colorScalePeriph, this.defaultParams.rPeriph,...
                this.defaultParams.gPeriph, this.defaultParams.bPeriph);

            this.shaderUniforms.colorScaleFov = glGetUniformLocation(this.shaderHandle3, 'colorScaleFov');
            glUniform3f(this.shaderUniforms.colorScaleFov, this.defaultParams.rFov,...
                this.defaultParams.gFov, this.defaultParams.bFov);

            this.shaderUniforms.slope = glGetUniformLocation(this.shaderHandle3, 'slope');
            glUniform1f(this.shaderUniforms.slope, this.defaultParams.slope);

            this.shaderUniforms.nxtoslope = glGetUniformLocation(this.shaderHandle3, 'nxtoslope');
            glUniform1f(this.shaderUniforms.nxtoslope, this.defaultParams.nxtoslope);

            this.shaderUniforms.flatwidth = glGetUniformLocation(this.shaderHandle3, 'flatwidth');
            glUniform1f(this.shaderUniforms.flatwidth, this.defaultParams.flatwidth);

            this.shaderUniforms.coeff = glGetUniformLocation(this.shaderHandle3, 'coeff');
            glUniform1f(this.shaderUniforms.coeff, this.defaultParams.coeff);

            this.shaderUniforms.invgamma = glGetUniformLocation(this.shaderHandle3, 'invgamma');
            glUniform1f(this.shaderUniforms.invgamma, this.defaultParams.invgamma);


            glUseProgram(0);           

            Screen('HookFunction', this.myOperator, 'AppendShader', 'UserDefinedBlit', 'MyShader3', this.shaderHandle3, []);
            
            if bitand(Screen('HookFunction', this.myOperator, 'ImagingMode'), mor(kPsychNeed16BPCFloat, kPsychNeed32BPCFloat)) == 0
                % Not yet set. Choose highest precision:
                Screen('HookFunction', this.myOperator, 'ImagingMode', mor(kPsychNeed32BPCFloat, Screen('HookFunction', this.myOperator, 'ImagingMode')));
                if debug > 3
                    fprintf('Add2DSeparableConvolutionToGLOperator: Increasing precision of operator to 32bpc float.\n');
                end
            end


        end

        function updateGLUniform(this,shaderHandle,uniformHandle,newUniformVal,dtype)

            glUseProgram(shaderHandle);

            switch dtype
                case "float"

                    switch length(newUniformVal)
                        case 1
                            glUniform1f(uniformHandle, newUniformVal);
                        case 2
                            glUniform2f(uniformHandle, newUniformVal(1),newUniformVal(2));
                        case 3
                            glUniform3f(uniformHandle, newUniformVal(1),newUniformVal(2),newUniformVal(3));
                        case 4
                            glUniform4f(uniformHandle, newUniformVal(1),newUniformVal(2),newUniformVal(3),newUniformVal(4));
                    end

                case "int"
                    glUniform1i(uniformHandle, newUniformVal);

            end
            
            glUseProgram(0);

        end

    end

end