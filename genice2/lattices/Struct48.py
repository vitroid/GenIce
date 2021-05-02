# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (24,44,8,0,)
"""
from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"SpaceFullerene": 'Sikiric 2010'},
        "usage": "No options available.",
        "brief": "A space fullerene."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.pairs = """
        120 392
        216 57
        124 330
        97 180
        92 191
        322 278
        97 5
        145 101
        253 228
        16 146
        284 215
        352 418
        229 35
        48 399
        176 92
        291 52
        378 109
        225 14
        333 55
        182 140
        12 18
        201 375
        91 50
        333 48
        311 238
        231 239
        108 401
        387 195
        289 363
        54 60
        25 169
        128 286
        100 350
        379 292
        11 291
        180 356
        264 418
        73 361
        174 282
        368 136
        2 403
        26 84
        13 71
        35 58
        92 120
        297 372
        43 304
        290 195
        387 263
        124 400
        56 209
        160 204
        169 185
        61 3
        422 370
        223 105
        166 170
        125 398
        271 228
        110 375
        158 119
        104 205
        49 320
        329 191
        174 146
        60 252
        364 425
        35 56
        89 265
        20 187
        118 52
        175 234
        67 32
        418 317
        426 120
        342 181
        95 341
        20 171
        368 380
        21 172
        100 262
        203 173
        134 217
        423 384
        226 126
        79 160
        360 345
        346 360
        86 191
        108 42
        403 426
        120 113
        377 15
        227 146
        196 63
        179 205
        269 115
        216 308
        251 385
        405 8
        219 322
        334 75
        18 159
        107 147
        64 132
        226 374
        292 249
        233 190
        308 116
        117 328
        117 329
        302 61
        316 312
        193 415
        372 140
        283 289
        341 335
        124 257
        21 175
        138 7
        408 115
        211 433
        391 190
        256 169
        360 384
        318 142
        127 300
        405 274
        225 314
        379 274
        386 99
        131 70
        1 40
        127 76
        196 203
        265 238
        354 258
        301 167
        95 165
        103 235
        256 67
        130 268
        47 183
        319 41
        106 38
        13 412
        274 346
        277 364
        365 279
        184 107
        50 10
        244 34
        101 160
        71 245
        237 59
        211 114
        71 266
        21 429
        51 376
        72 267
        112 306
        305 37
        55 77
        97 273
        137 49
        57 319
        51 306
        122 234
        321 238
        296 287
        334 215
        381 9
        119 375
        253 88
        200 138
        331 190
        182 391
        236 37
        237 267
        264 141
        168 166
        314 195
        220 105
        201 429
        394 311
        198 338
        334 381
        76 240
        61 344
        63 347
        197 260
        162 280
        227 289
        75 82
        322 326
        209 232
        229 298
        429 392
        119 43
        53 380
        131 90
        134 428
        118 18
        178 143
        24 151
        1 88
        197 397
        87 96
        162 259
        25 60
        298 248
        46 293
        50 173
        297 301
        371 412
        241 326
        381 115
        368 4
        358 125
        342 89
        230 427
        272 139
        373 398
        277 335
        97 399
        297 280
        172 375
        166 276
        231 182
        250 159
        13 189
        44 402
        45 422
        393 311
        200 45
        168 84
        19 210
        103 153
        270 147
        64 332
        36 114
        272 7
        173 430
        294 424
        17 30
        278 232
        356 406
        138 68
        243 316
        428 101
        357 431
        353 190
        220 407
        244 323
        104 58
        349 287
        207 147
        172 228
        433 41
        112 176
        168 361
        193 55
        30 395
        12 136
        341 371
        202 15
        144 396
        183 242
        150 382
        165 360
        420 427
        118 4
        112 38
        291 337
        69 48
        285 303
        434 176
        149 239
        184 214
        354 37
        261 431
        393 390
        301 210
        181 393
        81 239
        82 240
        78 8
        369 222
        143 5
        84 111
        304 72
        356 48
        367 208
        154 186
        201 88
        348 206
        374 170
        136 382
        383 106
        332 5
        399 77
        377 313
        358 168
        17 240
        98 223
        2 59
        22 146
        217 30
        285 252
        23 269
        314 166
        310 234
        155 204
        142 38
        154 78
        281 433
        57 102
        31 337
        377 401
        224 385
        319 361
        358 41
        56 68
        108 219
        216 73
        137 96
        433 398
        254 56
        237 306
        431 87
        367 157
        397 362
        284 233
        151 63
        133 365
        192 200
        8 292
        262 239
        38 392
        241 305
        226 362
        252 339
        121 424
        115 420
        148 397
        427 123
        353 85
        404 10
        286 66
        362 290
        327 422
        51 294
        397 169
        125 362
        340 384
        123 372
        157 16
        47 71
        355 311
        336 171
        400 364
        388 96
        244 9
        224 188
        403 303
        33 54
        421 7
        75 395
        2 33
        198 432
        258 161
        348 307
        416 411
        347 389
        313 42
        154 3
        152 401
        404 430
        404 432
        428 155
        218 256
        162 140
        201 113
        108 417
        423 412
        132 153
        226 276
        340 344
        386 177
        313 407
        328 54
        213 200
        407 22
        12 247
        33 163
        85 100
        318 214
        315 285
        267 234
        194 19
        295 321
        270 175
        15 402
        139 95
        250 26
        139 93
        217 171
        268 317
        216 53
        212 242
        32 185
        203 19
        416 243
        134 275
        117 315
        309 378
        349 336
        330 3
        265 247
        229 307
        423 425
        268 152
        379 83
        118 281
        155 240
        167 203
        145 287
        296 66
        99 365
        43 147
        27 293
        202 174
        5 128
        379 90
        315 263
        288 173
        122 74
        161 42
        269 151
        372 34
        197 111
        348 370
        405 186
        89 396
        373 116
        355 380
        342 321
        39 403
        386 81
        229 205
        236 402
        344 186
        309 406
        59 376
        103 393
        218 163
        39 121
        363 10
        11 181
        8 425
        302 421
        260 25
        227 44
        7 257
        106 267
        313 174
        329 121
        111 32
        383 411
        189 110
        352 232
        1 110
        231 233
        206 298
        31 114
        219 359
        424 163
        25 324
        94 45
        186 257
        273 222
        230 23
        410 294
        389 430
        242 271
        408 123
        28 351
        300 66
        14 84
        404 223
        135 67
        314 125
        262 82
        390 333
        409 282
        419 259
        236 326
        207 6
        284 85
        236 80
        46 320
        419 167
        14 148
        363 16
        354 198
        405 384
        65 293
        164 129
        149 255
        183 172
        91 196
        268 366
        86 122
        0 104
        91 65
        156 37
        135 252
        11 109
        434 304
        354 44
        406 396
        299 280
        213 421
        300 217
        310 191
        29 83
        371 400
        193 129
        80 202
        98 16
        254 213
        67 414
        408 62
        47 207
        235 337
        2 86
        351 87
        270 214
        76 233
        150 102
        4 26
        235 109
        331 30
        135 243
        13 40
        20 204
        21 207
        411 413
        343 26
        343 225
        6 266
        320 105
        112 410
        225 373
        377 156
        254 352
        351 299
        276 414
        99 160
        428 262
        52 380
        414 195
        133 246
        410 383
        150 109
        194 351
        390 143
        284 82
        99 336
        148 159
        409 258
        349 273
        11 355
        0 302
        322 152
        264 94
        28 151
        399 133
        250 374
        327 192
        134 350
        183 189
        158 142
        370 188
        29 107
        422 248
        387 324
        47 29
        366 282
        241 417
        130 348
        242 292
        130 305
        64 222
        367 283
        3 58
        312 324
        96 283
        211 116
        275 177
        127 331
        409 80
        126 312
        318 6
        339 263
        359 94
        312 290
        153 378
        62 140
        238 36
        344 138
        79 287
        391 9
        253 74
        101 81
        164 187
        220 80
        65 347
        353 395
        153 415
        320 338
        156 338
        251 370
        76 81
        23 334
        302 141
        90 107
        144 36
        106 426
        144 337
        149 75
        35 192
        27 105
        376 163
        126 260
        307 327
        237 434
        9 62
        170 73
        279 79
        4 361
        327 141
        220 156
        310 306
        132 325
        364 83
        427 215
        413 316
        184 266
        28 357
        342 356
        154 165
        224 366
        130 161
        350 255
        164 369
        340 199
        94 188
        218 316
        126 435
        329 33
        419 62
        214 88
        208 338
        374 398
        309 193
        318 429
        128 246
        124 95
        221 388
        149 244
        178 378
        264 213
        250 197
        179 141
        24 167
        317 278
        222 79
        279 129
        345 412
        50 261
        1 271
        192 61
        15 208
        232 359
        40 212
        265 382
        435 263
        328 413
        274 277
        199 68
        83 245
        142 72
        98 49
        68 139
        90 271
        261 388
        162 420
        70 346
        257 277
        221 430
        325 406
        315 416
        158 176
        350 395
        189 249
        288 347
        40 184
        218 60
        352 205
        391 255
        281 73
        177 255
        150 396
        199 165
        137 389
        27 10
        131 110
        373 159
        328 324
        180 178
        77 279
        119 6
        304 74
        0 340
        91 299
        394 53
        369 204
        407 258
        251 417
        256 290
        52 247
        209 298
        210 299
        388 363
        121 285
        235 394
        44 282
        85 420
        343 170
        63 194
        411 303
        78 400
        273 286
        93 104
        371 245
        187 365
        39 410
        325 55
        28 230
        31 281
        230 280
        39 92
        357 123
        330 93
        31 53
        241 206
        321 333
        331 177
        332 164
        136 319
        432 46
        295 325
        251 161
        243 435
        179 359
        296 246
        54 303
        70 249
        309 69
        117 51
        23 323
        74 113
        198 22
        137 194
        24 408
        65 87
        78 345
        426 122
        367 49
        181 178
        70 245
        89 291
        144 295
        353 381
        12 116
        27 157
        132 180
        209 45
        394 102
        86 434
        288 431
        289 432
        111 276
        199 58
        416 294
        179 248
        296 129
        42 366
        435 414
        417 248
        221 98
        272 335
        36 102
        275 336
        343 308
        409 385
        14 185
        358 148
        202 157
        64 77
        145 275
        402 152
        339 32
        349 133
        301 357
        423 335
        310 392
        231 34
        421 330
        196 259
        223 22
        0 272
        114 247
        418 307
        385 326
        206 278
        317 188
        215 34
        341 346
        389 46
        131 266
        369 66
        390 415
        24 288
        297 323
        376 413
        210 261
        212 29
        187 128
        43 228
        270 72
        300 155
        69 143
        387 185
        227 208
        103 295
        17 386
        283 293
        59 383
        260 339
        212 425
        20 17
        57 211
        345 249
        286 171
        69 246
        127 145
        368 308
        158 113
        305 401
        419 323
        182 100
        224 219
        221 19
        135 424
        254 93
        253 175
        18 41
        355 382
        269 259
        415 332
        """

        self.waters = """
        0.5 0.75 0.375
        0.69195 0.75 0.49227
        0.875 0.93807 0.59861
        0.625 0.25 0.38177
        0.8125 0.56695 0.76144
        0.625 0.56193 0.9014
        0.18305 0.4375 0.51144
        0.18807 0.93693 0.38149
        0.69195 0.06193 0.43405
        0.68305 0.125 0.04504
        0.0 0.94195 0.16456
        0.81307 0.56695 0.8288
        0.6875 0.05805 0.77707
        0.5 0.5625 0.4652
        0.6875 0.25 0.7152
        0.3125 0.55805 0.22294
        0.31193 0.94195 0.18405
        0.5 0.375 0.97712
        0.8125 0.93305 0.76144
        0.5 0.875 0.13177
        0.5 0.5 0.95544
        0.5 0.25 0.52288
        0.625 0.06695 0.20496
        0.5 0.625 0.06562
        0.75 0.25 0.09846
        0.68693 0.93305 0.67121
        0.6875 0.56695 0.73856
        0.0 0.75 0.17788
        0.375 0.43807 0.09861
        0.30805 0.0625 0.47294
        0.68305 0.4375 0.98856
        0.25 0.75 0.79456
        0.81193 0.44195 0.68405
        0.81193 0.06307 0.61851
        0.31695 0.125 0.04504
        0.81193 0.375 0.34861
        0.31193 0.05805 0.81595
        0.8125 0.43305 0.23856
        0.30805 0.56193 0.56595
        0.875 0.56193 0.59861
        0.5 0.75 0.47712
        0.0 0.05805 0.75774
        0.5 0.94195 0.25774
        0.0 0.125 0.52288
        0.125 0.25 0.22712
        0.0 0.875 0.32913
        0.81307 0.43305 0.17121
        0.375 0.25 0.48438
        0.68807 0.06307 0.88149
        0.5 0.625 0.17087
        0.0 0.93693 0.13812
        0.31307 0.25 0.61188
        0.75 0.75 0.79456
        0.375 0.56695 0.79504
        0.68693 0.93807 0.63149
        0.31193 0.06307 0.88149
        0.0 0.5 0.34846
        0.1875 0.25 0.7848
        0.68693 0.43807 0.36851
        0.125 0.93807 0.59861
        0.81193 0.875 0.6514
        0.68693 0.06193 0.36851
        0.80805 0.06193 0.06595
        0.68807 0.56307 0.11851
        0.25 0.75 0.90154
        0.0 0.56307 0.13812
        0.0 0.43693 0.94058
        0.0 0.375 0.67087
        0.0 0.625 0.36823
        0.81307 0.25 0.88812
        0.0 0.5 0.45544
        0.30805 0.4375 0.47294
        0.18305 0.875 0.54504
        0.1875 0.56695 0.76144
        0.81695 0.875 0.54504
        0.5 0.625 0.02288
        0.19195 0.25 0.99227
        0.375 0.93807 0.9014
        0.625 0.25 0.42087
        0.19195 0.93807 0.93405
        0.0 0.75 0.22712
        0.31695 0.0625 0.98856
        0.31695 0.5625 0.01144
        0.18305 0.125 0.45496
        0.8125 0.44195 0.72294
        0.0 0.5625 0.0348
        0.81695 0.06307 0.5788
        0.18693 0.43807 0.13149
        0.625 0.75 0.51563
        0.81307 0.93305 0.8288
        0.0 0.875 0.47712
        0.0 0.75 0.125
        0.81695 0.43693 0.5788
        0.31307 0.43807 0.36851
        0.18807 0.94195 0.31595
        0.125 0.43807 0.4014
        0.31193 0.375 0.1514
        0.75 0.75 0.90154
        0.5 0.875 0.17087
        0.5 0.0 0.95544
        0.0 0.75 0.02288
        0.19195 0.9375 0.97294
        0.18693 0.25 0.80942
        0.31193 0.625 0.84861
        0.5 0.56307 0.36188
        0.81307 0.75 0.19058
        0.375 0.75 0.57913
        0.18305 0.9375 0.48856
        0.5 0.75 0.27288
        0.0 0.44195 0.83544
        0.81695 0.5625 0.48856
        0.875 0.56695 0.70496
        0.18305 0.43693 0.5788
        0.81695 0.625 0.54504
        0.375 0.93305 0.79504
        0.80805 0.43807 0.06595
        0.5 0.125 0.76563
        0.5 0.25 0.625
        0.875 0.75 0.77288
        0.0 0.375 0.52288
        0.69195 0.56193 0.56595
        0.81193 0.43693 0.61851
        0.69195 0.93807 0.56595
        0.125 0.25 0.07913
        0.25 0.25 0.40154
        0.1875 0.05805 0.72294
        0.375 0.75 0.68439
        0.0 0.25 0.97712
        0.68305 0.43693 0.92121
        0.19195 0.25 0.91456
        0.625 0.25 0.27288
        0.0 0.625 0.47712
        0.125 0.75 0.88177
        0.68305 0.06307 0.92121
        0.875 0.75 0.98438
        0.0 0.5 0.65154
        0.8125 0.25 0.7848
        0.5 0.5 0.15154
        0.0 0.875 0.36823
        0.18807 0.56307 0.38149
        0.0 0.93693 0.05942
        0.5 0.94195 0.33544
        0.18305 0.625 0.54504
        0.68807 0.43693 0.88149
        0.18693 0.93305 0.8288
        0.0 0.0625 0.9652
        0.375 0.06695 0.20496
        0.18305 0.0625 0.51144
        0.8125 0.05805 0.72294
        0.5 0.875 0.02288
        0.0 0.25 0.82212
        0.625 0.43807 0.09861
        0.3125 0.43305 0.26144
        0.18693 0.56193 0.86851
        0.75 0.25 0.40154
        0.19195 0.5625 0.97294
        0.6875 0.55805 0.22294
        0.18693 0.75 0.19058
        0.0 0.5 0.54456
        0.6875 0.93305 0.73856
        0.31695 0.875 0.95496
        0.6875 0.06695 0.26144
        0.0 0.75 0.07212
        0.0 0.125 0.63177
        0.31695 0.43693 0.92121
        0.875 0.43807 0.4014
        0.1875 0.44195 0.72294
        0.625 0.06193 0.09861
        0.0 0.375 0.73438
        0.81193 0.05805 0.68405
        0.3125 0.56695 0.73856
        0.68305 0.625 0.95496
        0.69195 0.25 0.50774
        0.81307 0.06193 0.13149
        0.3125 0.94195 0.22294
        0.5 0.0625 0.5348
        0.0 0.375 0.56562
        0.68305 0.0625 0.98856
        0.81307 0.56193 0.86851
        0.5 0.75 0.32212
        0.875 0.75 0.88177
        0.68807 0.625 0.84861
        0.0 0.9375 0.0348
        0.625 0.25 0.48438
        0.30805 0.75 0.49227
        0.68693 0.25 0.69058
        0.875 0.06193 0.4014
        0.5 0.375 0.93439
        0.125 0.06695 0.29504
        0.69195 0.4375 0.47294
        0.875 0.25 0.01563
        0.69195 0.25 0.58544
        0.81193 0.125 0.34861
        0.18693 0.25 0.88812
        0.5 0.625 0.13177
        0.31307 0.25 0.69058
        0.81307 0.75 0.11188
        0.75 0.75 0.70544
        0.75 0.25 0.20544
        0.81193 0.56307 0.38149
        0.0 0.0 0.34846
        0.69195 0.5625 0.52707
        0.1875 0.75 0.2152
        0.68807 0.93693 0.11851
        0.31695 0.625 0.95496
        0.5 0.55805 0.33544
        0.875 0.43305 0.29504
        0.30805 0.25 0.50774
        0.375 0.43305 0.20496
        0.0 0.625 0.32913
        0.31193 0.93693 0.11851
        0.3125 0.05805 0.77707
        0.5 0.9375 0.4652
        0.18807 0.125 0.34861
        0.375 0.75 0.51563
        0.31695 0.375 0.04504
        0.3125 0.44195 0.77707
        0.80805 0.5625 0.97294
        0.0 0.0 0.65154
        0.3125 0.75 0.2848
        0.8125 0.75 0.2152
        0.5 0.0 0.15154
        0.125 0.75 0.92087
        0.68807 0.94195 0.18405
        0.1875 0.94195 0.27707
        0.5 0.25 0.72712
        0.25 0.75 0.70544
        0.25 0.25 0.20544
        0.81695 0.0625 0.51144
        0.68693 0.43305 0.3288
        0.31695 0.56307 0.0788
        0.19195 0.0625 0.02707
        0.18807 0.55805 0.31595
        0.125 0.25 0.01563
        0.5 0.06307 0.55942
        0.18693 0.56695 0.8288
        0.0 0.55805 0.24227
        0.18305 0.06307 0.5788
        0.5 0.125 0.82913
        0.31695 0.9375 0.01144
        0.31695 0.4375 0.98856
        0.8125 0.55805 0.27707
        0.69195 0.0625 0.47294
        0.18807 0.625 0.6514
        0.5 0.0 0.04456
        0.18305 0.375 0.45496
        0.80805 0.25 0.91456
        0.625 0.93305 0.79504
        0.68693 0.75 0.30942
        0.81695 0.375 0.45496
        0.625 0.75 0.72712
        0.8125 0.94195 0.27707
        0.81193 0.625 0.6514
        0.69195 0.9375 0.52707
        0.18807 0.375 0.34861
        0.68305 0.9375 0.01144
        0.0 0.125 0.67087
        0.125 0.06193 0.4014
        0.8125 0.06695 0.23856
        0.80805 0.75 0.08544
        0.625 0.75 0.68439
        0.18693 0.06193 0.13149
        0.19195 0.75 0.00774
        0.5 0.44195 0.66456
        0.31307 0.06695 0.3288
        0.68807 0.05805 0.81595
        0.18305 0.5625 0.48856
        0.30805 0.93807 0.56595
        0.375 0.25 0.27288
        0.68305 0.56307 0.0788
        0.30805 0.9375 0.52707
        0.81695 0.9375 0.48856
        0.31307 0.75 0.38812
        0.875 0.75 0.92087
        0.0 0.875 0.43439
        0.80805 0.9375 0.97294
        0.125 0.56695 0.70496
        0.18305 0.93693 0.42121
        0.125 0.43305 0.29504
        0.31695 0.06307 0.92121
        0.19195 0.75 0.08544
        0.125 0.75 0.77288
        0.1875 0.06695 0.23856
        0.18693 0.43305 0.17121
        0.19195 0.4375 0.02707
        0.68693 0.56193 0.63149
        0.80805 0.56193 0.93405
        0.0 0.06307 0.94058
        0.875 0.25 0.11823
        0.125 0.25 0.18439
        0.18807 0.05805 0.68405
        0.875 0.75 0.81562
        0.81695 0.125 0.45496
        0.0 0.55805 0.16456
        0.18807 0.43693 0.61851
        0.31193 0.875 0.84861
        0.0 0.25 0.92788
        0.31695 0.93693 0.0788
        0.81193 0.55805 0.31595
        0.18693 0.75 0.11188
        0.0 0.4375 0.9652
        0.375 0.06193 0.09861
        0.5 0.93693 0.36188
        0.625 0.75 0.61823
        0.0 0.0 0.54456
        0.6875 0.43305 0.26144
        0.30805 0.25 0.58544
        0.625 0.25 0.31562
        0.5 0.375 0.76563
        0.0 0.25 0.875
        0.5 0.25 0.57212
        0.5 0.375 0.82913
        0.31307 0.93305 0.67121
        0.5 0.875 0.23438
        0.3125 0.25 0.7152
        0.5 0.43693 0.63812
        0.18807 0.875 0.6514
        0.25 0.25 0.29456
        0.30805 0.5625 0.52707
        0.0 0.25 0.77288
        0.68807 0.55805 0.18405
        0.5 0.0 0.84846
        0.1875 0.55805 0.27707
        0.5 0.875 0.06562
        0.5 0.05805 0.66456
        0.18693 0.93807 0.86851
        0.0 0.625 0.26563
        0.68693 0.06695 0.3288
        0.5 0.06307 0.63812
        0.68693 0.25 0.61188
        0.375 0.25 0.38177
        0.80805 0.25 0.99227
        0.375 0.56193 0.9014
        0.5 0.125 0.86823
        0.5 0.5 0.04456
        0.30805 0.75 0.41456
        0.68305 0.875 0.95496
        0.125 0.75 0.81562
        0.625 0.43305 0.20496
        0.68693 0.56695 0.67121
        0.68693 0.75 0.38812
        0.18305 0.56307 0.42121
        0.68807 0.875 0.84861
        0.5 0.44195 0.74227
        0.81193 0.93693 0.38149
        0.69195 0.43807 0.43405
        0.0 0.625 0.43439
        0.81307 0.43807 0.13149
        0.75 0.25 0.29456
        0.80805 0.93807 0.93405
        0.80805 0.75 0.00774
        0.31193 0.56307 0.11851
        0.31307 0.43305 0.3288
        0.80805 0.4375 0.02707
        0.875 0.25 0.22712
        0.68807 0.44195 0.81595
        0.81307 0.93807 0.86851
        0.25 0.25 0.09846
        0.0 0.125 0.73438
        0.31307 0.75 0.30942
        0.81695 0.56307 0.42121
        0.0 0.44195 0.75774
        0.125 0.93305 0.70496
        0.18693 0.06695 0.17121
        0.30805 0.06193 0.43405
        0.5 0.125 0.93439
        0.3125 0.06695 0.26144
        0.31193 0.55805 0.18405
        0.6875 0.44195 0.77707
        0.19195 0.56193 0.93405
        0.875 0.06695 0.29504
        0.30805 0.43807 0.43405
        0.19195 0.06193 0.06595
        0.5 0.05805 0.74227
        0.375 0.75 0.72712
        0.81695 0.4375 0.51144
        0.18807 0.06307 0.61851
        0.5 0.625 0.23438
        0.0 0.43693 0.86188
        0.0 0.0 0.45544
        0.625 0.56695 0.79504
        0.68305 0.375 0.04504
        0.81307 0.25 0.80942
        0.25 0.75 0.59846
        0.69195 0.75 0.41456
        0.0 0.875 0.26563
        0.5 0.125 0.97712
        0.5 0.25 0.67788
        0.31193 0.125 0.1514
        0.68807 0.375 0.1514
        0.5 0.375 0.86823
        0.80805 0.0625 0.02707
        0.5 0.43693 0.55942
        0.5 0.5 0.84846
        0.31193 0.44195 0.81595
        0.68305 0.5625 0.01144
        0.0 0.05805 0.83544
        0.875 0.93305 0.70496
        0.3125 0.93305 0.73856
        0.625 0.93807 0.9014
        0.375 0.25 0.42087
        0.5 0.55805 0.25774
        0.1875 0.43305 0.23856
        0.75 0.75 0.59846
        0.81307 0.06695 0.17121
        0.81695 0.93693 0.42121
        0.0 0.06307 0.86188
        0.6875 0.94195 0.22294
        0.875 0.25 0.07913
        0.0 0.94195 0.24227
        0.125 0.56193 0.59861
        0.375 0.75 0.61823
        0.5 0.56307 0.44058
        0.31307 0.93807 0.63149
        0.18807 0.44195 0.68405
        0.31193 0.43693 0.88149
        0.31307 0.56193 0.63149
        0.6875 0.75 0.2848
        0.375 0.25 0.31562
        0.68305 0.93693 0.0788
        0.0 0.56307 0.05942
        0.31307 0.06193 0.36851
        0.81193 0.94195 0.31595
        0.5 0.75 0.42788
        0.0 0.375 0.63177
        0.5 0.93693 0.44058
        0.625 0.75 0.57913
        0.19195 0.43807 0.06595
        0.125 0.75 0.98438
        0.5 0.4375 0.5348
        0.68807 0.125 0.1514
        0.125 0.25 0.11823
        0.875 0.25 0.18439
        0.1875 0.93305 0.76144
        0.0 0.125 0.56562
        0.31307 0.56695 0.67121
        """

        self.coord = "relative"

        self.cages = """
        14 0.5 0.75 -0.17811
        15 0.5 0.75 -0.01482
        14 0.74773 0.25 1.65152
        14 0.5 0.99773 1.40152
        14 0.5 0.48219 0.29576
        14 -0.26781 -1.25 -0.95424
        12 0.0 -0.25 -0.15095
        15 0.5 -0.25 0.76482
        14 0.5 -0.25 0.92811
        15 0.0 -0.75 -0.51482
        12 0.5 0.5 0.5
        12 0.5 0.0 0.5
        14 -0.25227 -0.25 -0.65152
        15 0.0 0.75 0.51482
        12 1.0 0.5 1.0
        14 0.0 -0.01781 -0.20424
        12 -0.25 -1.25 -0.75
        12 0.5 0.25 0.34905
        12 0.0 -1.25 -0.70232
        14 0.0 0.25 0.42811
        14 0.0 1.00227 0.90152
        14 0.75227 0.75 1.15152
        12 0.5 0.75 0.09905
        14 0.5 0.49773 -0.40152
        14 0.76781 1.75 1.45424
        12 0.25 0.75 0.25
        12 0.25 1.25 0.75
        14 0.26781 1.25 0.95424
        14 0.25227 0.25 0.65152
        14 1.0 0.51781 0.79576
        14 0.0 0.75 0.625
        12 0.0 -0.75 -0.40095
        14 0.5 0.51781 0.70424
        14 1.0 1.49773 1.90152
        15 0.0 0.25 0.26482
        12 -0.25 -0.75 -0.25
        14 0.5 0.01781 0.29576
        14 0.5 0.25 0.125
        14 0.5 1.00227 -0.40152
        12 0.5 1.25 0.90095
        12 0.0 -0.25 0.95232
        14 0.24773 0.25 -0.15152
        14 0.0 -1.00227 -0.90152
        14 0.5 1.25 1.17811
        14 0.5 -0.25 0.875
        14 0.5 -0.01781 0.70424
        14 0.23219 1.25 0.54576
        12 -0.5 -0.25 -0.45232
        14 -0.23219 -0.75 -0.45424
        14 0.0 -0.75 -0.67811
        14 1.0 0.48219 1.20424
        14 -0.24773 0.25 -1.15152
        12 0.5 0.75 1.20232
        14 0.24773 -0.25 1.15152
        14 0.0 0.25 0.375
        12 0.0 0.25 0.04768
        14 0.26781 0.75 0.04576
        14 0.23219 0.75 0.45424
        14 0.5 1.50227 1.40152
        14 0.25227 -0.25 -0.65152
        12 0.5 1.25 -0.20232
        15 0.5 1.25 1.01482
        14 0.0 -0.25 -0.42811
        15 0.0 -0.25 -0.26482
        12 0.5 -0.25 0.65095
        14 0.0 0.01781 0.20424
        14 0.0 0.50227 -0.90152
        12 0.0 0.75 0.40095
        12 0.0 1.25 0.70232
        14 0.0 0.75 0.67811
        14 0.5 0.25 0.07189
        15 0.5 0.25 0.23518
        12 0.0 0.25 0.15095
        12 0.5 0.25 0.45232
        14 0.73219 1.25 0.95424
        12 0.0 0.0 0.0
        """

        self.bondlen = 3

        self.cell = """
        12.877303192396422 12.877303192396422 121.54450994600208
        """

        self.density = 0.6465960600253263

        self.cell = cellvectors(a=12.877303192396422,
                                b=12.877303192396422,
                                c=121.54450994600208)
