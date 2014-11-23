;(ql:quickload 'lispbuilder-sdl)
;(ql:quickload 'lispbuilder-sdl-gfx)

(defparameter *width* 100)
(defparameter *height* 100)
(defparameter *pixsize* 4) ;;One "pixel" in the sim is equal to 20um, this is an estimate, revise in stage two.
(defparameter *umperpix* 5)

(defparameter *cells* nil)
(defparameter *world* nil)
(defparameter *datafile* "Celldata.db")

(defun random-range (min max)
  (+ min (random (- max min) (make-random-state t))))

(defun new-cell (POS ATP NA AA FA G O2 CO2 DNA)
  (push (list :POS POS :ATP ATP :NA NA 
         :AA AA :FA FA :G G :O2 O2 :CO2 CO2 :DNA DNA) *cells*))

(defun new-world (TEMP PH NAC AAC FAC GC O2C CO2C RAD BPV)
  (setf *world* (list :TEMP TEMP :PH PH :NAC NAC 
                :AAC AAC :FAC FAC :GC GC :O2C O2 :CO2C CO2 :RAD RAD :BPV BPV)))

(defmacro fetch-value (accessor index lst)
  `(getf (nth ,index (reverse ,lst)) ,accessor))
  
(defun fetch-props (lst)
  (let ((nlst nil))
    (do ((i 0 (+ i 2)))
        ((>= i (length lst)))
      (push (nth i lst) nlst)) (reverse nlst)))  

(defun parse-gene (gene)
  (let ((gene-val 0))
    (dotimes (i (length gene))
      (setf gene-val (+ (* (expt 2 i) (nth i (reverse gene))) gene-val))) gene-val))
      
(defun parse-DNA (DNA-seq)
  (let ((DNA-val nil))
    (dotimes (i (length DNA-seq))
      (push (parse-gene (nth i DNA-seq)) DNA-val)) (reverse DNA-val)))

(defun dump-ci ()
  (dolist (entry *cells*)
    (format t "~%~{~a:~10t~a~%~}~%" entry)))

(defun save-ci (filename)
  (with-open-file (out filename
                   :direction :output
                   :if-exists :supersede)
    (with-standard-io-syntax
      (print *cells* out)
      (print *world* out))))

(defun load-ci (filename)
  (with-open-file (in filename)
    (with-standard-io-syntax
      (setf *cells* (read in))
      (setf *world* (read in)))))
      
(defun rand-mutate (DNA-seq)
  (let ((gene (random-range 0 (length DNA-seq))))
    (let ((base (random-range 0 (length (nth gene DNA-seq)))))
      (cond ((= (nth base (nth gene DNA-seq)) 0) (setf (nth base (nth gene DNA-seq)) 1))
            (t (setf (nth base (nth gene DNA-seq)) 0))) DNA-seq)))
            
(defun mutate-times (DNA-seq times)
  (dotimes (i times)
    (rand-mutate DNA-seq)))
  
(defun close-sim ()
  (format t "Saving simulation state and exiting...~%")
  (save-ci *datafile*)
  (exit))
  
(defun init-sim ()
  (format t "Loading simulation state and starting...~%")
  (when (probe-file *datafile*)
    (load-ci *datafile*))
  (format t "done."))

(defun Cell-Env (celli) ;;This is scientifically flawed, revise in stage two...
  (let ((prop nil) (cprop nil) (propcont nil) (cpropcont nil) (cellperinc nil))
    (do ((i 0 (1+ i)) (chance (random-range 0 100)))
        ((>= i 6)) ;;6 is for the four main macromolicules, make more portable later...
      
      (setf prop (nth i '(:NAC :AAC :FAC :GC :O2C :CO2C))) ;;make more portable later...
      (setf cprop (nth i '(:NA :AA :FA :G :O2 :CO2))) ;;make more portable later...
      (setf propcont (getf *world* prop))
      (setf cpropcont (fetch-value cprop celli *cells*))
      (setf cellperinc (+ (* (float (/ propcont (* 3 (* *width* *height*)))) (* *width* *height*)) cpropcont)) ;;Note: 3 is a placeholder for permeability
      
      (cond ((<= chance propcont) (setf (fetch-value cprop celli *cells*) cellperinc) (setf propcont (- propcont (/ cellperinc (* *width* *height*)))))
            (t (continue))))))
  
(defun Cell-Meta (celli)
  )  
  
(defun next-tick ()
  (dotimes (i (length *cells*))
  (Cell-Env i)
  (Cell-Meta)))

;(defun display-sim (&optional delay)
;  (sdl:with-init ()
;  (sdl:window (* *pixsize* *width*) (* *pixsize* *height*) :title-caption "Lisp Land")
;    (setf (sdl:frame-rate) 60)
;    (sdl:with-events ()
;      (:quit-event () t)
;      (:idle ()
;        (sleep delay)
;        (next-tick)
;        (dotimes (i (length *cells*))
;          (let ((POS (fetch-value :POS i *cells*)))
;                  (sdl-gfx:draw-box (sdl:rectangle :x (* (car POS) *pixsize*) :y (* (cdr POS) *pixsize*)
;                                    :w *pixsize* :h *pixsize*) :color sdl:*white*)))
;        (sdl:update-display)))))

(init-sim)
