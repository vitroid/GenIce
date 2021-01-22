# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (40,20,8,12,)
"""
desc={"ref": {"SpaceFullerene": 'Sikiric 2010'},
      "usage": "No options available.",
      "brief": "A space fullerene."
      }

import genice2.lattices
from genice2.cell import cellvectors

class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.pairs="""
        2 323
        223 93
        3 308
        429 439
        399 16
        385 338
        342 277
        24 419
        7 307
        24 280
        326 277
        409 160
        205 240
        64 455
        242 129
        74 62
        76 147
        233 189
        276 398
        451 348
        377 250
        125 408
        187 23
        265 348
        419 31
        313 88
        188 196
        54 216
        61 298
        192 49
        349 383
        26 333
        195 37
        225 442
        410 318
        255 93
        35 427
        334 145
        435 189
        230 221
        15 66
        247 408
        183 177
        165 196
        86 109
        36 384
        28 343
        389 196
        192 105
        40 77
        159 317
        61 141
        210 285
        368 404
        358 149
        304 329
        43 52
        87 396
        327 135
        120 39
        99 375
        325 443
        428 208
        377 379
        44 112
        43 130
        414 176
        316 288
        97 437
        152 393
        420 21
        306 426
        57 94
        174 223
        358 194
        357 193
        219 109
        10 14
        416 344
        383 184
        24 432
        421 216
        104 446
        67 415
        210 187
        145 330
        412 184
        111 13
        399 208
        149 148
        204 205
        391 224
        439 329
        14 278
        256 38
        301 317
        57 69
        307 269
        375 114
        266 28
        246 195
        59 273
        242 351
        74 362
        311 54
        132 422
        356 411
        445 30
        259 315
        446 30
        32 445
        394 333
        126 454
        282 17
        333 173
        1 325
        376 158
        121 306
        390 370
        359 186
        72 381
        132 278
        46 20
        257 151
        152 266
        373 159
        324 384
        126 212
        225 291
        103 156
        301 161
        152 185
        364 389
        122 153
        10 361
        404 300
        224 95
        258 18
        400 277
        22 12
        273 153
        199 210
        185 445
        201 355
        79 418
        80 419
        80 420
        57 394
        213 238
        246 432
        314 21
        417 79
        364 96
        383 3
        215 312
        390 9
        243 338
        386 177
        226 115
        259 192
        260 193
        221 163
        387 205
        139 371
        171 234
        56 12
        131 343
        453 366
        213 111
        377 170
        230 310
        345 6
        451 82
        266 211
        143 335
        266 430
        190 237
        0 49
        257 122
        1 51
        64 369
        156 115
        331 341
        349 279
        434 101
        281 424
        0 366
        402 75
        392 302
        393 305
        421 426
        248 440
        229 194
        67 214
        372 22
        306 90
        429 84
        392 128
        294 273
        279 447
        104 170
        399 436
        19 217
        27 118
        120 98
        339 144
        208 347
        173 318
        402 347
        261 320
        170 16
        42 249
        336 127
        257 362
        120 25
        164 97
        268 328
        88 439
        166 84
        334 434
        103 264
        174 5
        150 175
        255 151
        387 220
        234 73
        146 318
        247 42
        219 155
        412 168
        199 248
        201 249
        337 397
        281 403
        123 264
        15 360
        322 300
        6 65
        181 304
        117 68
        267 418
        241 425
        275 61
        80 140
        444 171
        378 52
        124 284
        315 188
        349 426
        106 233
        27 4
        354 32
        46 91
        395 198
        220 420
        399 384
        231 100
        141 13
        225 441
        194 109
        46 380
        271 29
        136 113
        136 115
        222 161
        280 127
        304 252
        313 133
        286 175
        406 337
        368 431
        421 292
        413 407
        94 268
        438 18
        177 8
        36 288
        207 100
        67 329
        373 425
        229 293
        382 113
        357 292
        414 355
        423 263
        340 98
        97 328
        154 344
        2 296
        34 147
        444 268
        404 4
        358 85
        381 397
        253 290
        367 424
        90 183
        373 141
        74 134
        69 40
        413 89
        142 300
        253 83
        158 272
        292 382
        96 406
        191 365
        429 374
        353 273
        372 450
        225 22
        350 141
        197 55
        39 41
        444 363
        150 407
        203 343
        274 15
        217 78
        218 290
        60 271
        351 228
        256 35
        362 215
        20 4
        107 338
        34 374
        63 232
        148 58
        231 265
        370 151
        91 427
        117 211
        200 274
        10 156
        129 45
        65 153
        198 220
        164 168
        199 424
        289 50
        358 99
        186 9
        454 33
        42 434
        137 186
        191 397
        376 347
        209 49
        27 100
        293 446
        401 352
        438 31
        155 133
        0 315
        114 75
        79 262
        193 426
        289 407
        138 9
        136 403
        112 380
        2 410
        334 312
        105 110
        207 380
        423 29
        353 338
        289 382
        219 350
        373 297
        276 232
        316 452
        408 251
        339 182
        68 17
        112 35
        385 425
        252 65
        402 72
        437 90
        314 303
        120 455
        103 8
        154 215
        81 375
        66 361
        200 54
        36 347
        117 131
        412 328
        104 208
        430 298
        272 81
        169 179
        411 189
        414 330
        327 45
        377 0
        409 116
        201 135
        54 279
        349 14
        241 166
        352 217
        309 370
        345 135
        29 439
        423 228
        261 162
        63 22
        286 70
        350 354
        388 302
        130 83
        379 406
        180 411
        293 32
        165 436
        53 355
        357 422
        125 102
        113 8
        56 450
        336 64
        129 255
        313 85
        31 182
        361 279
        60 155
        342 310
        430 13
        210 73
        102 237
        369 335
        331 107
        368 51
        369 12
        57 363
        393 297
        89 409
        367 156
        192 96
        409 70
        174 319
        41 232
        357 164
        339 113
        62 448
        119 449
        365 211
        94 440
        180 9
        2 324
        198 48
        25 235
        47 145
        19 38
        39 336
        451 78
        28 445
        128 133
        366 95
        56 24
        142 237
        34 29
        356 237
        238 203
        228 344
        291 340
        403 168
        239 364
        318 308
        421 183
        243 6
        417 85
        442 160
        59 252
        242 257
        299 250
        300 251
        367 328
        44 212
        108 317
        391 117
        214 181
        354 254
        7 290
        159 305
        314 89
        361 264
        149 81
        386 332
        378 12
        392 211
        446 58
        126 451
        247 87
        454 443
        223 145
        387 218
        41 450
        164 73
        55 325
        52 18
        431 101
        204 455
        376 178
        353 169
        394 187
        316 324
        322 326
        281 103
        173 389
        53 345
        321 431
        259 389
        437 285
        86 99
        434 233
        287 212
        221 371
        11 296
        106 197
        335 38
        43 303
        138 167
        249 93
        119 387
        280 303
        385 317
        261 323
        351 215
        146 196
        19 265
        129 416
        81 178
        198 127
        128 68
        124 310
        209 58
        11 453
        253 235
        199 168
        19 454
        227 92
        258 303
        367 440
        166 341
        124 44
        213 365
        418 214
        413 157
        436 453
        46 231
        55 33
        11 364
        185 114
        415 416
        138 139
        105 379
        307 256
        272 17
        91 205
        342 163
        368 44
        160 21
        26 268
        271 59
        378 143
        154 255
        241 107
        61 169
        121 323
        311 398
        308 77
        296 69
        187 299
        333 172
        418 295
        362 309
        428 343
        386 216
        438 274
        144 37
        283 16
        212 48
        150 63
        353 108
        97 171
        296 173
        452 234
        321 71
        91 401
        403 292
        189 309
        25 175
        158 379
        236 428
        363 299
        408 5
        429 302
        96 381
        174 42
        119 253
        136 123
        371 100
        221 167
        71 202
        236 95
        246 66
        452 162
        254 28
        249 312
        126 310
        282 149
        14 440
        165 239
        393 298
        90 422
        372 332
        344 107
        223 370
        165 283
        236 72
        217 98
        224 381
        116 182
        157 83
        202 284
        246 31
        243 134
        311 332
        287 51
        267 111
        430 295
        209 245
        280 420
        272 293
        144 123
        282 417
        315 363
        64 82
        395 218
        146 40
        200 50
        171 320
        407 70
        67 84
        226 447
        327 330
        441 235
        397 203
        271 181
        291 369
        40 323
        376 245
        36 239
        351 176
        48 82
        283 288
        340 235
        437 424
        60 108
        206 16
        276 18
        86 185
        132 3
        43 269
        188 69
        222 99
        448 190
        441 157
        286 442
        360 50
        4 163
        150 441
        1 207
        20 240
        62 294
        452 410
        147 228
        76 415
        45 176
        25 449
        209 170
        20 348
        431 284
        34 214
        258 37
        105 245
        241 415
        152 213
        119 21
        116 195
        372 286
        78 455
        306 278
        206 324
        444 172
        47 102
        23 110
        47 101
        258 432
        154 134
        226 383
        391 191
        360 332
        53 294
        201 87
        130 7
        435 390
        178 75
        275 108
        327 93
        433 176
        47 5
        453 188
        322 356
        337 131
        121 285
        139 251
        385 65
        7 427
        341 374
        197 124
        356 167
        92 135
        248 184
        448 227
        72 239
        400 265
        261 206
        378 232
        405 263
        197 277
        337 17
        242 76
        259 172
        282 191
        161 85
        398 450
        240 78
        106 319
        316 110
        92 396
        180 106
        386 70
        227 102
        137 190
        260 73
        354 375
        270 359
        26 184
        438 144
        388 417
        289 216
        321 5
        287 443
        195 264
        62 122
        238 114
        118 240
        137 411
        359 134
        186 396
        299 320
        55 71
        110 250
        206 250
        111 297
        207 251
        394 248
        104 158
        263 416
        162 283
        243 294
        231 352
        398 432
        405 252
        254 13
        71 319
        442 449
        218 118
        435 167
        245 288
        325 112
        155 262
        56 336
        275 219
        159 222
        331 147
        301 88
        262 161
        346 360
        448 355
        181 262
        26 308
        193 3
        94 77
        284 233
        254 238
        276 200
        274 447
        52 419
        400 371
        220 427
        230 202
        1 321
        447 123
        172 23
        388 267
        143 39
        224 49
        449 204
        422 346
        234 23
        435 202
        340 401
        412 115
        400 33
        148 391
        244 10
        388 365
        341 298
        433 59
        75 30
        244 132
        82 38
        413 339
        68 229
        285 320
        109 32
        137 74
        138 125
        128 79
        130 80
        313 302
        89 37
        146 162
        404 380
        222 350
        77 278
        395 287
        295 374
        281 183
        346 177
        270 122
        236 436
        334 309
        260 121
        35 48
        127 269
        405 345
        27 51
        227 330
        180 326
        311 66
        41 175
        433 423
        392 295
        142 101
        335 98
        194 133
        359 151
        270 6
        405 45
        88 60
        204 291
        319 390
        53 433
        139 326
        179 425
        230 33
        346 382
        263 329
        428 30
        86 305
        331 169
        322 163
        95 131
        178 58
        443 256
        260 410
        116 8
        148 229
        247 142
        395 307
        396 125
        267 84
        314 83
        87 190
        118 352
        270 92
        342 348
        414 312
        275 305
        11 384
        366 406
        157 140
        402 203
        244 15
        160 140
        179 153
        301 304
        166 297
        290 401
        143 269
        182 140
        244 226
        63 50
        76 179
        """

        self.waters="""
        0.19092 0.32156 0.9918
        0.69092 0.82156 0.9918
        0.5 0.25711 0.52095
        0.81592 0.17888 0.61187
        0.30908 0.83709 0.23528
        0.68708 0.7617 0.98993
        0.80908 0.6321 0.30857
        0.81592 0.93216 0.00536
        0.125 0.09622 0.15209
        0.81292 0.7383 0.48993
        0.18408 0.13514 0.79938
        0.5 0.30464 0.68044
        0.0 0.98574 0.61647
        0.5 0.48574 0.61647
        0.30908 0.16291 0.73528
        0.0 0.097 0.67553
        0.18792 0.31435 0.36116
        0.5 0.41417 0.01992
        0.69337 0.03331 0.80223
        0.75 0.88072 0.55608
        0.18408 0.86486 0.29938
        0.30663 0.98381 0.09917
        0.0 0.01246 0.49147
        0.68792 0.26156 0.11493
        0.18163 0.001 0.80301
        0.5 0.97141 0.38289
        0.81292 0.21181 0.83077
        0.5 0.84916 0.17458
        0.31837 0.43353 0.60457
        0.19337 0.57207 0.97879
        0.125 0.39714 0.46127
        0.0 0.04611 0.89708
        0.31837 0.44056 0.3498
        0.68792 0.83171 0.6186
        0.31837 0.55945 0.8498
        0.0 0.88914 0.89811
        0.68708 0.33127 0.48867
        0.5 0.05317 0.99468
        0.875 0.90378 0.65209
        0.5 0.96418 0.6035
        0.18792 0.23844 0.61493
        0.5 0.00023 0.55483
        0.0 0.73253 0.92793
        0.69337 0.97758 0.88333
        0.18408 0.83849 0.86932
        0.625 0.63945 0.95768
        0.0 0.87021 0.20611
        0.5 0.74289 0.02095
        0.18408 0.89724 0.80721
        0.0 0.34288 0.02195
        0.81592 0.06785 0.50536
        0.5 0.84288 0.02195
        0.81837 0.001 0.80301
        0.125 0.63745 0.15183
        0.5 0.09417 0.59946
        0.81292 0.81575 0.7452
        0.18163 0.98976 0.67801
        0.31208 0.24129 0.86493
        0.0 0.3825 0.20729
        0.375 0.59622 0.15209
        0.19337 0.54355 0.182
        0.31837 0.53116 0.49113
        0.31592 0.66152 0.36932
        0.81837 0.03116 0.49113
        0.18163 0.93353 0.60457
        0.68408 0.60276 0.30721
        0.18408 0.08489 0.73902
        0.68163 0.55945 0.8498
        0.31592 0.43216 0.00536
        0.31292 0.25868 0.73993
        0.375 0.05101 0.30666
        0.68792 0.78869 0.77478
        0.875 0.36255 0.65183
        0.68792 0.21131 0.27478
        0.19092 0.67844 0.4918
        0.875 0.39714 0.46127
        0.5 0.597 0.67553
        0.125 0.21038 0.67475
        0.375 0.89714 0.46127
        0.5 0.48755 0.99147
        0.0 0.98755 0.99147
        0.68408 0.41511 0.23902
        0.125 0.90378 0.65209
        0.69337 0.98381 0.09917
        0.80663 0.53331 0.80223
        0.80663 0.48381 0.09917
        0.0 0.47141 0.38289
        0.0 0.71274 0.1808
        0.0 0.53582 0.1035
        0.5 0.03582 0.1035
        0.18408 0.16152 0.36932
        0.0 0.903 0.17553
        0.68708 0.68426 0.2452
        0.81208 0.68565 0.86116
        0.18708 0.21181 0.83077
        0.19092 0.3679 0.80857
        0.68408 0.33849 0.86932
        0.0 0.19536 0.18044
        0.69337 0.92793 0.47879
        0.80663 0.46669 0.30223
        0.69092 0.83709 0.23528
        0.31292 0.7617 0.98993
        0.5 0.7327 0.17793
        0.25 0.11928 0.05608
        0.31592 0.36486 0.29938
        0.68408 0.32113 0.11187
        0.0 0.76731 0.67793
        0.0 0.59417 0.59946
        0.125 0.55101 0.30666
        0.19337 0.46669 0.30223
        0.625 0.28962 0.17475
        0.68163 0.48976 0.67801
        0.0 0.85012 0.92869
        0.875 0.09622 0.15209
        0.80663 0.42793 0.47879
        0.875 0.13945 0.95768
        0.18163 0.06647 0.10457
        0.19092 0.41953 0.88435
        0.5 0.8825 0.20729
        0.375 0.95923 0.17689
        0.5 0.94684 0.49468
        0.375 0.2098 0.42968
        0.5 0.64988 0.42869
        0.625 0.10286 0.96127
        0.18708 0.81575 0.7452
        0.68708 0.74132 0.23993
        0.375 0.86255 0.65183
        0.375 0.94899 0.80666
        0.31837 0.46884 0.99113
        0.68408 0.63514 0.79938
        0.81837 0.96884 0.99113
        0.31592 0.39724 0.80721
        0.0 0.15765 0.58109
        0.19337 0.48381 0.09917
        0.0 0.65712 0.52195
        0.81208 0.66829 0.1186
        0.75 0.11928 0.05608
        0.125 0.7098 0.42968
        0.68792 0.75871 0.36493
        0.81292 0.78819 0.33077
        0.0 0.01426 0.11647
        0.5 0.51246 0.49147
        0.18792 0.76156 0.11493
        0.69337 0.95645 0.682
        0.69337 0.07207 0.97879
        0.5 0.71332 0.92587
        0.0 0.25734 0.58502
        0.31592 0.58489 0.73902
        0.0 0.40584 0.09946
        0.80908 0.42656 0.12958
        0.69337 0.02242 0.38333
        0.68408 0.67888 0.61187
        0.0 0.46418 0.6035
        0.5 0.61086 0.39811
        0.0 0.65084 0.67458
        0.31837 0.51025 0.17801
        0.125 0.13945 0.95768
        0.81837 0.01025 0.17801
        0.5 0.37021 0.20611
        0.80663 0.52242 0.38333
        0.18163 0.01025 0.17801
        0.68163 0.51025 0.17801
        0.0 0.26748 0.42793
        0.31208 0.81435 0.36116
        0.81292 0.18426 0.2452
        0.0 0.31328 0.58618
        0.875 0.54077 0.67689
        0.5 0.76748 0.42793
        0.68792 0.16829 0.1186
        0.31592 0.56785 0.50536
        0.19092 0.33709 0.23528
        0.0 0.2327 0.17793
        0.81292 0.2617 0.98993
        0.68708 0.25868 0.73993
        0.81208 0.74129 0.86493
        0.5 0.00344 0.42983
        0.375 0.63945 0.95768
        0.18408 0.10276 0.30721
        0.80908 0.38904 0.30217
        0.5 0.58583 0.51992
        0.0 0.75711 0.52095
        0.5 0.54913 0.05624
        0.0 0.04913 0.05624
        0.30908 0.1321 0.30857
        0.68792 0.18565 0.86116
        0.0 0.44684 0.49468
        0.875 0.7098 0.42968
        0.5 0.24266 0.08502
        0.18792 0.28869 0.77478
        0.31208 0.73844 0.61493
        0.18792 0.71131 0.27478
        0.80908 0.41953 0.88435
        0.80908 0.32156 0.9918
        0.69092 0.17844 0.4918
        0.125 0.45923 0.17689
        0.30663 0.07207 0.97879
        0.0 0.28727 0.6808
        0.0 0.80464 0.68044
        0.30908 0.91953 0.88435
        0.5 0.18672 0.08618
        0.69092 0.07344 0.62958
        0.0 0.68672 0.08618
        0.5 0.78727 0.6808
        0.625 0.40378 0.65209
        0.18163 0.94056 0.3498
        0.18408 0.91511 0.23902
        0.31292 0.28819 0.33077
        0.81592 0.82113 0.11187
        0.375 0.36055 0.45768
        0.0 0.34916 0.17458
        0.5 0.21274 0.1808
        0.125 0.44899 0.80666
        0.30908 0.8679 0.80857
        0.80663 0.45645 0.682
        0.5 0.54611 0.89708
        0.19092 0.66291 0.73528
        0.5 0.09818 0.44507
        0.625 0.89714 0.46127
        0.5 0.90584 0.09946
        0.31837 0.499 0.30301
        0.18408 0.93216 0.00536
        0.5 0.80368 0.42633
        0.68163 0.499 0.30301
        0.68708 0.71181 0.83077
        0.0 0.37221 0.90027
        0.0 0.99177 0.36647
        0.81592 0.13514 0.79938
        0.5 0.69536 0.18044
        0.19092 0.61096 0.80217
        0.19092 0.42656 0.12958
        0.5 0.81328 0.58618
        0.81592 0.86486 0.29938
        0.69337 0.01619 0.59917
        0.18708 0.75868 0.73993
        0.81292 0.24132 0.23993
        0.69337 0.96669 0.30223
        0.125 0.36255 0.65183
        0.31292 0.74132 0.23993
        0.68163 0.43353 0.60457
        0.81208 0.33171 0.6186
        0.30908 0.88904 0.30217
        0.80908 0.57344 0.62958
        0.5 0.62979 0.70611
        0.0 0.62779 0.40027
        0.0 0.12979 0.70611
        0.80908 0.33709 0.23528
        0.18163 0.05945 0.8498
        0.0 0.74266 0.08502
        0.5 0.19632 0.92633
        0.0 0.69632 0.92633
        0.375 0.28962 0.17475
        0.875 0.78962 0.17475
        0.625 0.59622 0.15209
        0.625 0.95923 0.17689
        0.5 0.45087 0.55624
        0.80908 0.66291 0.73528
        0.81592 0.89724 0.80721
        0.5 0.65765 0.58109
        0.5 0.02859 0.88289
        0.875 0.2902 0.92968
        0.625 0.2098 0.42968
        0.18792 0.25871 0.36493
        0.5 0.51426 0.11647
        0.875 0.60286 0.96127
        0.375 0.10286 0.96127
        0.875 0.86055 0.45768
        0.19337 0.45645 0.682
        0.68163 0.501 0.80301
        0.0 0.21332 0.92587
        0.625 0.94899 0.80666
        0.68408 0.66152 0.36932
        0.31837 0.56647 0.10457
        0.5 0.403 0.17553
        0.31592 0.60276 0.30721
        0.81592 0.08489 0.73902
        0.19337 0.52242 0.38333
        0.625 0.04077 0.67689
        0.0 0.81388 0.52129
        0.18408 0.17888 0.61187
        0.5 0.1175 0.70729
        0.30663 0.97758 0.88333
        0.375 0.13745 0.15183
        0.68408 0.43216 0.00536
        0.0 0.30368 0.42633
        0.31208 0.78869 0.77478
        0.31208 0.21131 0.27478
        0.30663 0.02242 0.38333
        0.5 0.87221 0.90027
        0.81208 0.31435 0.36116
        0.69092 0.08047 0.38435
        0.69092 0.92656 0.12958
        0.0 0.95389 0.39708
        0.69092 0.1321 0.30857
        0.31592 0.41511 0.23902
        0.19092 0.6321 0.30857
        0.31837 0.501 0.80301
        0.5 0.26731 0.67793
        0.80663 0.51619 0.59917
        0.19337 0.51619 0.59917
        0.31208 0.26156 0.11493
        0.125 0.78962 0.17475
        0.80663 0.54355 0.182
        0.0 0.49656 0.92983
        0.5 0.99656 0.92983
        0.68163 0.56647 0.10457
        0.0 0.50344 0.42983
        0.30908 0.17844 0.4918
        0.69092 0.91953 0.88435
        0.875 0.21038 0.67475
        0.375 0.71038 0.67475
        0.31208 0.83171 0.6186
        0.30908 0.07344 0.62958
        0.18792 0.68565 0.86116
        0.0 0.49977 0.05483
        0.5 0.99977 0.05483
        0.125 0.2902 0.92968
        0.68708 0.28819 0.33077
        0.875 0.55101 0.30666
        0.81208 0.23844 0.61493
        0.81292 0.75868 0.73993
        0.18708 0.24132 0.23993
        0.625 0.7902 0.92968
        0.18708 0.78819 0.33077
        0.31292 0.2383 0.48993
        0.5 0.28668 0.42587
        0.81592 0.83849 0.86932
        0.0 0.78668 0.42587
        0.68708 0.66873 0.98867
        0.0 0.18612 0.02129
        0.80663 0.57207 0.97879
        0.5 0.68612 0.02129
        0.19092 0.57344 0.62958
        0.18408 0.06785 0.50536
        0.68792 0.24129 0.86493
        0.31292 0.71181 0.83077
        0.81837 0.93353 0.60457
        0.30663 0.95645 0.682
        0.5 0.38914 0.89811
        0.0 0.59818 0.44507
        0.81837 0.06647 0.10457
        0.81837 0.94056 0.3498
        0.125 0.54077 0.67689
        0.18708 0.83127 0.48867
        0.375 0.40378 0.65209
        0.0 0.6175 0.70729
        0.875 0.63745 0.15183
        0.0 0.11086 0.39811
        0.625 0.36055 0.45768
        0.125 0.86055 0.45768
        0.5 0.15084 0.67458
        0.5 0.49177 0.36647
        0.31592 0.63514 0.79938
        0.69092 0.88904 0.30217
        0.19092 0.58047 0.38435
        0.5 0.45389 0.39708
        0.18792 0.66829 0.1186
        0.31208 0.75871 0.36493
        0.81592 0.16152 0.36932
        0.875 0.45923 0.17689
        0.80908 0.67844 0.4918
        0.0 0.08583 0.51992
        0.30908 0.11096 0.80217
        0.31592 0.67888 0.61187
        0.18708 0.2617 0.98993
        0.68708 0.31575 0.7452
        0.875 0.44899 0.80666
        0.31592 0.33849 0.86932
        0.18708 0.16873 0.98867
        0.30908 0.82156 0.9918
        0.0 0.95087 0.55624
        0.625 0.71038 0.67475
        0.68792 0.81435 0.36116
        0.18163 0.03116 0.49113
        0.68163 0.53116 0.49113
        0.19337 0.53331 0.80223
        0.68163 0.44056 0.3498
        0.68408 0.36486 0.29938
        0.31592 0.32113 0.11187
        0.81837 0.98976 0.67801
        0.5 0.34236 0.08109
        0.0 0.84236 0.08109
        0.80908 0.3679 0.80857
        0.81592 0.10276 0.30721
        0.69092 0.16291 0.73528
        0.5 0.31388 0.52129
        0.80908 0.58047 0.38435
        0.30908 0.08047 0.38435
        0.30908 0.92656 0.12958
        0.80663 0.47758 0.88333
        0.81208 0.28869 0.77478
        0.68792 0.73844 0.61493
        0.0 0.40183 0.94507
        0.19337 0.47758 0.88333
        0.0 0.50023 0.55483
        0.5 0.23253 0.92793
        0.5 0.90183 0.94507
        0.81208 0.71131 0.27478
        0.68408 0.39724 0.80721
        0.375 0.04077 0.67689
        0.31292 0.33127 0.48867
        0.81292 0.83127 0.48867
        0.81592 0.91511 0.23902
        0.75 0.38072 0.55608
        0.625 0.13745 0.15183
        0.18408 0.82113 0.11187
        0.75 0.61928 0.05608
        0.5 0.35012 0.92869
        0.625 0.05101 0.30666
        0.81208 0.76156 0.11493
        0.30663 0.04355 0.182
        0.68708 0.2383 0.48993
        0.18708 0.7383 0.48993
        0.81292 0.16873 0.98867
        0.69337 0.04355 0.182
        0.31292 0.66873 0.98867
        0.68408 0.58489 0.73902
        0.80908 0.61096 0.80217
        0.68163 0.46884 0.99113
        0.5 0.50823 0.86647
        0.0 0.00823 0.86647
        0.18163 0.96884 0.99113
        0.5 0.12779 0.40027
        0.0 0.14988 0.42869
        0.125 0.60286 0.96127
        0.31208 0.16829 0.1186
        0.68408 0.56785 0.50536
        0.5 0.15712 0.52195
        0.0 0.91417 0.01992
        0.25 0.38072 0.55608
        0.0 0.52859 0.88289
        0.31837 0.48976 0.67801
        0.375 0.7902 0.92968
        0.30663 0.03331 0.80223
        0.25 0.61928 0.05608
        0.18792 0.74129 0.86493
        0.5 0.75734 0.58502
        0.18792 0.33171 0.6186
        0.18708 0.18426 0.2452
        0.81837 0.05945 0.8498
        0.0 0.55317 0.99468
        0.31208 0.18565 0.86116
        0.81837 0.999 0.30301
        0.18163 0.999 0.30301
        0.69092 0.8679 0.80857
        0.0 0.24289 0.02095
        0.19337 0.42793 0.47879
        0.19092 0.38904 0.30217
        0.69092 0.11096 0.80217
        0.31292 0.68426 0.2452
        0.30663 0.96669 0.30223
        0.30663 0.01619 0.59917
        0.25 0.88072 0.55608
        0.81208 0.25871 0.36493
        0.31292 0.31575 0.7452
        0.625 0.86255 0.65183
        0.30663 0.92793 0.47879
        """

        self.coord= "relative"

        self.cages="""
        15 0.0 0.40814 0.68695
        12 0.27348 -0.47751 0.99363
        12 0.22652 0.02249 0.99363
        12 0.22652 0.97751 1.49363
        14 1.0 1.14149 0.18373
        12 0.0 0.32921 0.81306
        14 0.5 0.4019 0.43005
        14 -0.26369 0.37574 1.05053
        12 0.0 -0.42418 0.80451
        16 0.5 0.55343 0.30763
        12 0.27348 0.47751 0.49363
        14 0.26369 -0.37574 0.55053
        16 0.5 0.20918 0.67366
        15 0.5 0.09186 0.18695
        15 0.5 0.45863 1.17099
        14 0.23631 0.12426 0.55053
        14 0.26369 0.37574 0.05053
        16 0.0 -0.29082 0.67366
        12 0.25167 0.28478 0.55404
        12 0.24833 0.21522 0.05404
        15 0.0 0.04137 0.67099
        14 -0.5 -0.14566 1.42358
        16 0.0 0.29082 0.17366
        12 0.0 0.49357 0.73978
        16 0.0 -0.20952 0.92996
        14 0.23631 0.87574 1.05053
        12 0.0 0.17254 0.79334
        14 -0.26369 -0.37574 0.55053
        12 0.5 0.92418 0.30451
        12 1.0 0.50643 0.23978
        12 -0.27348 -0.47751 0.99363
        16 0.5 0.44657 0.80763
        12 0.0 0.25029 0.80205
        12 0.25167 0.71522 1.05404
        12 -0.5 -0.00643 0.73978
        12 0.74833 0.28478 0.55404
        14 1.0 0.64566 0.92358
        16 0.5 0.79082 0.17366
        12 1.0 0.67079 0.31306
        14 0.5 0.14566 0.92358
        15 -0.5 -0.09186 0.68695
        14 0.5 0.64149 0.18373
        12 0.5 0.24971 0.30205
        12 0.77348 0.02249 0.99363
        16 0.0 -0.05343 0.80763
        16 0.5 0.70952 0.42996
        12 -0.5 -0.24971 0.80205
        14 1.0 0.9019 0.43005
        12 0.5 0.07582 0.80451
        14 0.76369 0.12426 0.55053
        14 -0.5 -0.4019 0.93005
        12 1.0 0.55232 0.49228
        12 0.5 0.00643 0.23978
        14 0.0 -0.14149 0.68373
        14 0.76369 0.87574 0.05053
        12 0.5 0.05232 0.49228
        15 0.5 0.54137 0.67099
        12 -0.27348 0.47751 1.49363
        14 0.0 0.35434 1.42358
        12 0.77348 0.97751 0.49363
        16 1.0 1.20952 0.42996
        12 -0.5 -0.17079 0.81306
        12 0.74833 0.71522 0.05404
        12 1.0 0.82746 0.29334
        12 1.0 0.74971 0.30205
        16 0.0 0.05343 0.30763
        14 0.0 0.0981 0.93005
        12 -0.5 -0.32746 0.79334
        12 0.5 0.32746 0.29334
        12 0.0 0.42418 0.30451
        14 -0.5 -0.64149 0.68373
        15 1.0 0.59186 0.18695
        15 0.0 -0.04137 1.17099
        12 -0.24833 0.21522 1.05404
        16 -0.5 -0.70952 0.92996
        12 0.24833 -0.21522 0.55404
        12 0.0 0.44768 0.99228
        12 -0.24833 -0.21522 0.55404
        12 0.5 0.17079 0.31306
        12 -0.5 -0.05232 0.99228
        """

        self.bondlen = 3


        self.cell = """
        15.106429533183366 98.03113793959956 21.936018274703592
        """

        self.density = 0.41957820843744786



        self.cell = cellvectors(a=15.106429533183366,
                           b=98.03113793959956,
                           c=21.936018274703592)
