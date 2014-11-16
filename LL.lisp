;(ql:quickload 'lispbuilder-sdl)
;(ql:quickload 'lispbuilder-sdl-gfx)

(defparameter *width* 100)
(defparameter *height* 100)
(defparameter *pixsize* 5)

(defparameter *cells* nil)
(defparameter *world* nil)
(defparameter *datafile* "Celldata.db")

(defun random-range (min max)
  (+ min (random (- max min) (make-random-state t))))

(defun new-cell (POS ATP NA AA FA G DNA)
  (push (list :POS POS :ATP ATP :NA NA 
         :AA AA :FA FA :G G :DNA DNA) *cells*))

(defun new-world (TEMP PH NAC AAC FAC GC RAD BPV)
  (setf *world* (list :TEMP TEMP :PH PH :NAC NAC 
                :AAC AAC :FAC FAC :GC GC :RAD RAD :BPV BPV)))

(defmacro fetch-value (accessor index lst)
  `(getf (nth ,index (reverse ,lst)) ,accessor))

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

(defun next-tick ()
  ;;(Cel-Env) Function for reacting to the enviroment, for now, absorbing macromolicules.
  )

;(defun display-sim ()
;  (sdl:with-init ()
;  (sdl:window (* *pixsize* *width*) (* *pixsize* *height*) :title-caption "Lisp Land")
;    (setf (sdl:frame-rate) 60)
;    (sdl:with-events ()
;      (:quit-event () t)
;      (:idle ()
;        ;;(next-tick)
;        (dotimes (i (length *cells*))
;          (let ((POS (fetch-value :POS i *cells*)))
;                  (sdl-gfx:draw-box (sdl:rectangle :x (* (car POS) *pixsize*) :y (* (cdr POS) *pixsize*)
;                                    :w *pixsize* :h *pixsize*) :color sdl:*white*)))
;        (sdl:update-display)))))

(init-sim)
