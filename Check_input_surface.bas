Attribute VB_Name = "Check_input_surface"
Option Explicit



Public Sub LayThongTin()


    '---------------------------------Lay thong tin tinh sut bo---------------------------------
    '//______________________________________________________________________________________//
    '//______________________________________________________________________________________//
    '//______________________________________________________________________________________//

    Public c As Integer, phi As Integer, u  As Integer
    c = 5000 'kN
    phi = 20 'do
    u = 0 'N
    Public omega As Single
    omega = 0
    Public Delta As Integer
    Delta = 10
    Public D As Single
    D = 0
    Public Number_of_Center As Integer
    Number_of_Center = 1
    Public x_center_entry As Integer, y_center_entry As Integer
    x_center_entry = 90 'm
    y_center_entry = 55 'm
    Public RL_entry As Integer
    RL_entry = 1 'm
    Public cr As Integer
    cr = Round(Math.Sqr(Namber_of_Center))
    'Mat cat xet
    Public h() As Single
    ReDim h(15)
    h(1) = 45
    h(2) = 45
    h(3) = 42.5
    h(4) = 40
    h(5) = 37.5
    h(6) = 35
    h(7) = 32.5
    h(8) = 30
    h(9) = 27.5
    h(10) = 25
    h(11) = 22.5
    h(12) = 20
    h(13) = 17.5
    h(14) = 15
    h(15) = 15




    Public rl() As Single
    Public entry_point As Single, exit_point As Single
    Public R As Single



    Public alpha() As Single
    Public beta() As Single
    Public R As Single
    Public x() As Single
    Public f As Single
    Public kc_d As Single
    Public FS_m() As Single
    Public FS_f() As Single


    Public lamda() As Single




    Public Matrix_C() As Integer
    Public N() As Single
    Public W() As Single

    Public Number_of_Slices As Integer
    'Bien chay
    Public i As Integer, j As Integer


    Public Tolerance as Single
    Tolerance = 0.1


End Sub


