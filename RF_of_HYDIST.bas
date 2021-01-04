Attribute VB_Name = "FS"


Option Explicit
Type Surface
    first As Integer
    final As Integer
    c As Single
    alpha(first To final) As Single
    beta(first To final) As Single
    phi As Single
    u As Single
    omega As Single
    x(first To final) As Single
    aa(first To final) As Single
    W(first To final) As Single
    fx(first to final) as Single
    kW As Single
    D As Single
    A As Single  
End Type
 
 Type Normal_Force
    N_f() as Single
    N_m() as Single
End Type




'============================================================================================================
'============================================================================================================


Public Sub Defining_Center()
    'Tao luoi tam
    ReDim Matrix_C(Number_of_Center, 2)
    If Number_of_Center = 1 Then`
        Matrix_C(1, 1) = x_center_entry
        Matrix_C(1, 2) = y_center_entry
    Else
        For j = 1 To cr + 1
            For i = 1 To cr
                Dim dem As Integer
                dem = (j - 1) * cr + i
                If dem <= Number_of_Center Then
                    Matrix_C(dem, 1) = (i - 1) * Delta / 10 + x_center_entry
                    Matrix_C(dem, 2) = (j - 1) * Delta / 10 + y_center_entry
                End If
            Next
        Next
    End If
End Sub

'------------------------------------------------------------------------------------s-
Public Sub Defining_Radius_Lines()
    'Tao Duong xac dinh ban kinh
    'Thuong se nam trong mien dat
    'Radius Lines chi dung cac duong nam ngang ( horrizontal)
    '----------------------------------
    ReDim rl(cr + 1)
    For i = 1 To cr + 1
        rl(i) = RL_entry + (i - 1) * Delta / 100
    Next
    
End Sub

'-----------------------------------------------------------------------------------------
Public Sub Defining_Entry_Exit_Point_Radius(ByVal rl As Single, ByVal c1 As Integer, ByVal c2 As Integer, ByRef entry_point As Single, ByRef exit_point As Single, ByRef R As Single)
    R = Math.Abs(c2 - rl)
    i = 1
    Dim tam_l As Single, tam_r As Single
    Do
        i = i + 1
        tam_l = ((i - 1) * 10 - c1) ^ 2 + (h(i) - c2) ^ 2
        tam_r = (i * 10 - c1) ^ 2 + (h(i) - c2) ^ 2
    Loop Until (tam_l < R ^ 2 And tam_r > R ^ 2)
    entry_point = Math.Sqr(R ^ 2 - (h(i) - c2) ^ 2) + c1
    If entry_point < (i - 1) * 10 Or entry_point > i * 10 Then entry_point = c1 - Math.Sqr(R ^ 2 - (h(i) - c2) ^ 2)
    
    Do
        i = i + 1
        tam_l = ((i - 1) * 10 - c1) ^ 2 + (h(i) - c2) ^ 2
        tam_r = (i * 10 - c1) ^ 2 + (h(i) - c2) ^ 2
    Loop Until (tam_l < R ^ 2 And tam_r > R ^ 2)
    exit_point = Math.Sqr(R ^ 2 - (h(i) - c2) ^ 2) + c1
    If exit_point < (i - 1) * 10 Or entry_point > i * 10 Then exit_point = c1 - Math.Sqr(R ^ 2 - (h(i) - c2) ^ 2)

End Sub









'---------------------------------------------------------------------------------------------
'Chia lai slice
Public Sub Descretizating_Slices(ByVal rl As Single, ByVal c1 As Single, ByVal c2 As Single)

    'If UBound(h) > 30 Then Number_of_Slices = UBound(h)
    'Else
            'Noi suy slice
    'End If
    
    Number_of_Slices = UBound(h)
    Call Defining_Entry_Exit_Point_Radius(rl, c1, c2, entry_point, exit_point, R)
    
    
    
    
    
End Sub

'===========================================================================================================================================================
'===========================================================================================================================================================
'===========================================================================================================================================================
'===========================================================================================================================================================

Private Sub Calculating_N(ByVal XrXl As Boolean, ByVal F_m As Boolean, ByVal FS_f_m As Single,ByVal lamda As Single, ByVal MC as Surface,byRef N as Normal_Force )
    Dim count as Integer
    Dim Xl_Xr(MC.first to MC.final) as Single
    ReDim N.N_f(MC.first to MC.final) as Single
    ReDim N.N_m(MC.first to MC.final) as Single
    For count = MC.first To MC.final
        If XrXl Then
            Dim tam1 As Single, tam2 As Single
            tam1 = (MC.c * MC.beta(count) - MC.u * MC.beta(count)) / FS_f_m
            tam2 = N(count) * (Math.Tan(MC.phi) * Math.Cos(MC.alpha(count)) / FS_f_m - Math.Sin(MC.alpha(count)))
            Xr_Xl(count) = MC.fx(count) * lamda * (tam1 + tam2 - MC.kW + MC.D * Math.Cos(MC.omega))
        Else
            Xr_Xl(count) = 0
        End If
        Dim tam3 As Single
        Dim tam4 As Single
        tam3 = (MC.c * MC.beta(count) * Math.Sin(MC.alpha(count)) + MC.u * MC.beta(count) * Math.Tan(MC.phi)) / FS_f_m
        tam4 = Math.Cos(MC.alpha(count)) + Math.Sin(MC.alpha(count)) * Math.Tan(MC.phi) / FS_f_m
        If F_m Then
            N.N_f(count) = (MC.W(count) + Xr_Xl(count) - tam3) / tam4
        Else
            N.N_m(count) = (MC.W(count) + Xr_Xl(count) - tam3) / tam4
        End If

    Next
        
End Sub

'_____________________________________________________________________________________________________________________________________________________________________________





Private Function Calculating_FS(ByVal N as Normal_Force, ByVal F_m As Boolean, ByVal MC as Surface) As Single
        Dim count As Integer
        Dim tam1 As Single, tam2 As Single, tam3 As Single
        Dim FS As Single
        tam1 = 0
        tam2 = 0
        tam3 = 0
        If F_m Then
            For count = MC.first To MC.final
                tam1 = tam1 + MC.c * MC.beta(count) * Math.Cos(MC.alpha(count))
                tam2 = tam2 + (N.N_f(count) - MC.u * MC.beta(count)) * Math.Tan(MC.phi) * Math.Cos(MC.alpha(count))
                tam3 = tam3 + N.N_f(count) * Math.Sin(MC.alpha(count)) - MC.A
            Next
            FS = (tam1 + tam2) / tam3
        Else
            For count = MC.first To MC.final
                tam1 = tam1 + MC.c * MC.beta(count) * R
                tam2 = tam2 + (N.N_m(count) - MC.u * MC.beta(count)) * R * Math.Tan(MC.phi)
                tam3 = tam3 + MC.W(count) * MC.x(count) - MC.A * MC.aa(count)
            Next
            FS = (tam1 + tam2) / tam3
        End If
        Calculating_FS = FS
End Function




'/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Public Sub FoS_Calculating(ByRef fx() As Single, ByRef lamda() As Single)
    
    
    
    '////////////////////////////////////////////////////
    '///////////////////////////////////////////////////
    ' MAIN CALCULATING BELOW
    ' Building an plot with every surface in vertical of lamda changing series and A list of value F_f,F_m
    'After got that, we take intersection point of two lines F_f(lamda) and F_m(lamda).
    'At that piont, we got F_m = F_f = F.O.S of the current surface
    '   With every lamda values, We use RAPID SOLVER - Multi Iterations to have Convergence Status
    '   The result of RAPID SOLVER is list of F_m corresponding with lamda and also with F_f
    '


    'MAIN 4 STAGES of FOS CALCULATING
    'Main Concept of RAPID SOLVER - Newton Raphson Algorithm
    '-----------------------------------------------------------------
    '   STAGE 1: 1 Iteration : Set X = E = 0 at every slice ; set lamda = 0 ;
    '   STAGE 2: 4-6 Iterations :    X = E = 0 as XlXr = False in Calculating_N ; to convergence F_m and F_f
    '   STAGE 3: Main RAPID SOLVER - According to Newton Raphson Algorithm
    '       set initial lamda = 2/3 slope between crest and toe
    '       Selected Tolerance
    '       Iterations until we got Abs(F_f - F_m) < Tolerance
    '   STAGE 4: We got a list of F_f(lamda) and F_m(lamda)
    '       Aferthat, we plot F_f(lamda) and F_m(lamda) ; the intersection is value of F.O.S that we finding




    Dim iii As Integer
    Dim jjj As Integer
    Dim kkk As Integer

    For kkk = 1 To UBound(lamda)
        Dim lamda_crr As Single
        lamda_crr = lamda(kkk)
        Dim FS_f_crr As Single
        Dim FS_m_crr As Single

        For iii = 1 To UBound(fx)

            Dim fx_curr As Single

            fx_curr = fx(iii)

        'STAGE 1:
        '//////////////////////////////////////////////////

            FS_f_crr = Calculating_FS(True, first, final)
            FS_m_crr = Calculating_FS(False, first, final)
        
        'STAGE 2
        '////////////////////////////////
        
            For jjj = 1 To 6
                Call Calculating_N(False, True, FS_f_crr, first, final, fx_curr, lamda_crr)
                FS_f_crr = Calculating_FS(True, first, final)
            Next
            For jjj = 1 To 6
                Call Calculating_N(False, False, FS_m_crr, first, final, fx_curr, lamda_crr)
                FS_m_crr = Calculating_FS(False, first, final)
            Next
        
        
        'STAGE 3
        '/////////////////////////////////
            Do
                Call Calculating_N(True, True, FS_f_crr, first, final, fx_curr, lamda_crr)
                FS_f_crr = Calculating_FS(True, first, final)
                Call Calculating_N(True, False, FS_m_crr, first, final, fx_curr, lamda_crr)
                FS_m_crr = Calculating_FS(False, first, final)
            Loop While (Math.Abs(FS_f_crr - FS_m_crr) > Tolerance)
       Next
       FS_f(kkk) = FS_f_crr
       FS_m(kkk) = FS_m_crr

        
    Next



    '////////////////////////////////////////////////////////////////////////////
    'End fx,lamda loop
    '
    'STAGE 4: F_f(lamda) and F_m(lamda)
    
    Dim note As Integer
    For kkk = 1 To (UBound(lamda) - 1)
        If (FS_f(kkk + 1) - FS_f(kkk)) * (FS_m(kkk + 1) - FS_m(kkk)) < 0 Then
            note = kkk
            Exit For
        End If
    Next
    'Interpolating FOS
    Dim tam As Single
    Dim tam1 As Single
    Dim lamda_final As Single
    tam = (FS_f(note + 1) - FS_m(note + 1) - FS_f(note) + FS_m(note)) / (lamda(note + 1) - lamda(note))
    tam1 = lamda(note) * (FS_f(note + 1) - FS_f(note) - FS_m(note + 1) + FS_m(note)) / (lamda(note + 1) - lamda(note)) + FS_m(note) - FS_f(note)
    lamda_final = tam1 / tam
    FoS_Calculating = lamda_final * (FS_f(note + 1) - FS_f(note)) / (lamda(note + 1) - lamda(note)) - lamda(note) * (FS_f(note + 1) - FS_f(note)) / (lamda(note + 1) - lamda(note)) + FS_f(note)

End Sub
